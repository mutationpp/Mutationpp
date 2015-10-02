/**
 * @file CollisionDB.cpp
 *
 * @brief Implementation of CollisionDB type.
 */

/*
 * Copyright 2014-2015 von Karman Institute for Fluid Dynamics (VKI)
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

#include "CollisionDBNew.h"
#include "Species.h"
#include "Thermodynamics.h"
#include "XMLite.h"
#include "Utilities.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;

#include <cassert>
#include <fstream>
#include <iostream>
using namespace std;

namespace Mutation {
    namespace Transport {

//==============================================================================

CollisionDBNew::CollisionDBNew(
    const string& db_name, const Thermodynamics::Thermodynamics& thermo) :
    m_database(databaseFileName(db_name, "transport")),
    m_thermo(thermo),
    m_tabulate(true)
{
    // Loop over the species and create the list of species pairs
    const vector<Species>& species = m_thermo.species();
    XmlElement& root = m_database.root();

    // Determine if we should tabulate collision integrals if possible
    root.getAttribute("tabulate", m_tabulate, m_tabulate);
    cout << "tabulating: " << m_tabulate << endl;

    for (int i = 0; i < species.size(); ++i)
        for (int j = i; j < species.size(); ++j)
            m_pairs.push_back(CollisionPairNew(species[i], species[j], &root));
}

//==============================================================================

const CollisionGroup& CollisionDBNew::group(const string& name)
{
    // Minimal check on the string argument
    assert(isValidGroupName(name));

    // Check if this group is already being managed
    map<string, CollisionGroup>::iterator iter = m_groups.find(name);
    if (iter != m_groups.end())
        return iter->second.update(
            (name[3] == 'e' ? m_thermo.Te() : m_thermo.T()), m_thermo);

    // Create a new group to manage this type
    CollisionGroup& new_group = m_groups.insert(
        make_pair(name, CollisionGroup(m_tabulate))).first->second;

    string type  = name.substr(0,3); // Q11, Q22, Bst, ...
    string pairs = name.substr(3);   // ee, ei, ii, ij

    const int ns = m_thermo.nSpecies();
    const int e  = (m_thermo.hasElectrons() ? 1 : 0);
    const int k  = e*ns;

    std::vector<CollisionPairNew>::iterator start, end;
    std::vector<CollisionPairNew> diag;

    if (pairs == "ee") {
        // just electron-electron pair
        start = m_pairs.begin();
        end   = m_pairs.begin()+e;
    } else if (pairs == "ei") {
        // electron-heavy pairs
        start = m_pairs.begin();
        end   = m_pairs.begin()+k;
    } else if (pairs == "ij") {
        // all heavy-heavy pairs
        start = m_pairs.begin()+k;
        end   = m_pairs.end();
    } else if (pairs == "ii") {
        // diagonal heavy-heavy pairs
        // create a temporary list of the heavy diagonal components
        for (int i = 0, index = k; i < ns-e; index += ns-e-i, i++)
            diag.push_back(m_pairs[index]);
        start = diag.begin();
        end   = diag.end();
    }

    new_group.manage(start, end, &CollisionPairNew::get, type);

    // Just call the function again (don't recode temperature selection, etc.)
    return (*this).group(name);
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation

