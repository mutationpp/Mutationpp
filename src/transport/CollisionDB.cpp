/**
 * @file CollisionDB.cpp
 *
 * @brief Implementation of CollisionDB type.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "CollisionDB.h"
#include "Species.h"
#include "Thermodynamics.h"
#include "XMLite.h"
#include "Utilities.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;

#include <Eigen/Dense>
using namespace Eigen;

#include <cassert>
#include <fstream>
#include <iostream>
using namespace std;

namespace Mutation {
    namespace Transport {

//==============================================================================

CollisionDB::CollisionDB(
    const string& db_name, const Thermodynamics::Thermodynamics& thermo) :
    m_database(databaseFileName(db_name, "transport")),
    m_thermo(thermo),
    m_ng(thermo.nGas()),
    m_nh(thermo.nHeavy() - thermo.nCondensed()),
    m_tabulate(false), m_table_min(300.),  m_table_max(20000.), m_table_del(100.),
    m_mass(m_ng),
    m_etai(m_nh), m_etafac(m_nh),
    m_nDei(m_ng*(thermo.hasElectrons() ? 1 : 0)),
    m_Deifac(m_ng*(thermo.hasElectrons() ? 1 : 0)),
    m_nDij(m_nh*(m_nh+1)/2),
    m_Dijfac(m_nh*(m_nh+1)/2),
    m_Dim(m_ng),
    m_L01ei(m_ng*(thermo.hasElectrons() ? 1 : 0)),
    m_L02ei(m_ng*(thermo.hasElectrons() ? 1 : 0))
{
    XmlElement& root = m_database.root();

    // Load global options
    XmlElement::const_iterator opts = root.findTag("global-options");
    if (opts != root.end()) {
        // Determine if we should tabulate collision integrals when possible
        XmlElement::const_iterator tabulate = opts->findTag("tabulate");
        if (tabulate != opts->end()) {
            m_tabulate = true;
            tabulate->getAttribute("Tmin", m_table_min, m_table_min);
            tabulate->getAttribute("Tmax", m_table_max, m_table_max);
            tabulate->getAttribute("dT",   m_table_del, m_table_del);
        }

        // Check the table data
        if (m_tabulate) {
            tabulate->parseCheck(m_table_min > 0.0, "Tmin must be positive.");
            tabulate->parseCheck(m_table_max > 0.0, "Tmax must be positive.");
            tabulate->parseCheck(m_table_del > 0.0, "dT must be positive.");
            tabulate->parseCheck(m_table_min < m_table_max, "Tmin must be > Tmax.");

            double size = (m_table_max-m_table_min)/m_table_del;
            tabulate->parseCheck(std::abs(size-int(size))/size < 1.0e-15,
                "(Tmax - Tmin)/dT must be a positive whole number.");
        }
    }

    // Loop over the species and create the list of species pairs
    const vector<Species>& species = m_thermo.species();
    for (int i = 0; i < nSpecies(); ++i)
        for (int j = i; j < nSpecies(); ++j)
            m_pairs.push_back(CollisionPair(species[i], species[j], &root));

    // Compute species masses
    for (int i = 0; i < nSpecies(); ++i)
        m_mass(i) = thermo.speciesMw(i) / NA;

    // Compute eta factors
    const int nh = nHeavy();
    const int k  = nSpecies()-nh;
    m_etafac = 5./16.*std::sqrt(PI*KB)*m_mass.tail(nh).sqrt();

    // Compute the nDei factors
    if (nSpecies() > nh) {
        m_Deifac.fill(3./16.*std::sqrt(TWOPI*KB/m_mass(0)));
        m_Deifac(0) *= 2./SQRT2;
    }

    // Compute the nDij factors
    for (int i = k, index = 0; i < nSpecies(); ++i)
        for (int j = i; j < nSpecies(); ++j, index++)
            m_Dijfac(index) = 3./16.*std::sqrt(TWOPI*KB*
                (m_mass(i)+m_mass(j))/(m_mass(i)*m_mass(j)));
}

//==============================================================================

int CollisionDB::nSpecies() const { return m_ng; }

//==============================================================================

int CollisionDB::nHeavy() const { return m_nh; }

//==============================================================================

Map<const ArrayXd> CollisionDB::X() const {
    return Eigen::Map<const ArrayXd>(thermo().X(), nSpecies());
}

//==============================================================================

Map<const ArrayXd> CollisionDB::Y() const {
    return Eigen::Map<const ArrayXd>(thermo().Y(), nSpecies());
}

//==============================================================================

CollisionDB::GroupType CollisionDB::groupType(const string& name)
{
    const int len = name.length();
    if (name[len-2] == 'e') {
        if (name [len-1] == 'e') return EE;
        if (name [len-1] == 'i') return EI;
    }
    else if (name[len-2] == 'i') {
        if (name [len-1] == 'i') return II;
        if (name [len-1] == 'j') return IJ;
    }
    return BAD_TYPE;
}

//==============================================================================


const CollisionGroup& CollisionDB::group(const string& name)
{
    // Minimal check on the string argument
    GroupType type = groupType(name);
    assert(type != BAD_TYPE);

    // Check if this group is already being managed
    map<string, CollisionGroup>::iterator iter = m_groups.find(name);
    if (iter != m_groups.end()) return iter->second.update(
        (type < II ? m_thermo.Te() : m_thermo.T()), m_thermo);

    // Create a new group to manage this type
    CollisionGroup& new_group = m_groups.insert(
        make_pair(name, CollisionGroup(
            m_tabulate, m_table_min, m_table_max, m_table_del))).first->second;

    // Determine start and end iterator for this group
    const int ns = nSpecies();
    const int e  = (m_thermo.hasElectrons() ? 1 : 0);
    const int k  = e*ns;

    std::vector<CollisionPair>::iterator start = m_pairs.begin();
    std::vector<CollisionPair>::iterator end   = m_pairs.end();
    std::vector<CollisionPair> diag;

    switch (type) {
    case EE: end = start + e; break;
    case EI: end = start + k; break;
    case IJ: start += k; break;
    case II:
        // create a temporary list of the heavy diagonal components
        for (int i = 0, index = k; i < ns-e; index += ns-e-i, i++)
            diag.push_back(m_pairs[index]);
        start = diag.begin();
        end   = diag.end();
        break;
    default:
        cout << "Bad collision integral group type: '"
             << name.substr(name.length()-2) << "' in group name: '"
             << name << "'." << endl;
        cout << "Allowed group types are 'ee', 'ei', 'ii', and 'ij'." << endl;
        exit(1);
    }

    // Manage the collision integrals
    string kind = name.substr(0, name.length()-2);
    new_group.manage(start, end, &CollisionPair::get, kind);

    // Compute integrals and return the group
    return new_group.update(
        (type < II ? m_thermo.Te() : m_thermo.T()), m_thermo);
}

//==============================================================================

const ArrayXd& CollisionDB::etai()
{
    return (m_etai = std::sqrt(m_thermo.T()) * m_etafac / Q22ii());
}

//==============================================================================

const ArrayXd& CollisionDB::nDei()
{
    if (m_nDei.size() > 0)
        m_nDei = std::sqrt(m_thermo.Te()) * m_Deifac / Q11ei();
    return m_nDei;
}

//==============================================================================

const ArrayXd& CollisionDB::nDij()
{
    return (m_nDij = std::sqrt(m_thermo.T()) * m_Dijfac / Q11ij());
}

//==============================================================================

const ArrayXd& CollisionDB::Dim(bool include_one_minus_x)
{
    const int ns = nSpecies();
    const int nh = m_nh;
    const int k  = ns - nh;
    const ArrayXd X = Map<const ArrayXd>(m_thermo.X(), ns).max(1.0e-12);

    // Electron
    if (k > 0) {
        nDei(); // update nDei
        m_Dim(0) = (X / m_nDei).tail(nh).sum();
        m_Dim.tail(nh) = X(0) / m_nDei.tail(nh);
    } else
        m_Dim.setZero();

    // Heavies
    nDij(); // update nDij
    for (int i = 0, index = 1; i < nh; ++i, ++index) {
        for (int j = i+1; j < nh; ++j, ++index) {
            m_Dim(i+k) += X(j+k)/m_nDij(index);
            m_Dim(j+k) += X(i+k)/m_nDij(index);
        }
    }

    if (include_one_minus_x) {
        // Compute (1-X_i) as sum_j!=i X_j for better accuracy
        for (int i = 0; i < ns; ++i)
            m_Dim(i) = (X.head(i).sum() + X.tail(ns-i-1).sum()) / m_Dim(i);
    }
    else {
        m_Dim = 1.0 / m_Dim;
    }

    // Remove number density and return
    return (m_Dim /= m_thermo.numberDensity());
}

//==============================================================================

const ArrayXd& CollisionDB::L01ei()
{
    if (m_L01ei.size() > 0) {
        const double theta = m_thermo.Te() / m_thermo.T();
        m_L01ei = theta*X()*(3.*Q12ei()-2.5*Q11ei());
        m_L01ei(0) = -m_L01ei.tail(nHeavy()).sum() / theta;
    }
    return m_L01ei;
}

//==============================================================================

const ArrayXd& CollisionDB::L02ei()
{
    if (m_L02ei.size() > 0) {
        const double theta = m_thermo.Te() / m_thermo.T();
        m_L02ei = theta*X()*(10.5*Q12ei()-6.*Q13ei()-35./8.*Q11ei());
        m_L02ei(0) = -m_L02ei.tail(nHeavy()).sum() / theta;
    }
    return m_L02ei;
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation

