/**
 * @file MillikanWhite.cpp
 *
 * @brief Implementation of classes related to Millikan and White model.
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


#include "MillikanWhite.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Mutation::Thermodynamics;

namespace Mutation {
    namespace Transfer {


struct MillikanWhiteModelData::Impl
{
    size_t index = 0;
    double mw = 0.0;
    double omegav = 1.0E-20; // m^2
    Eigen::ArrayXd a;
    Eigen::ArrayXd b;

    Impl(size_t size) : a(size), b(size) { }
};


MillikanWhiteModelData::MillikanWhiteModelData(
    const class Thermodynamics& thermo, size_t index, double thetav
) :
    m_impl(std::make_unique<Impl>(thermo.nHeavy()))
{
    m_impl->index = index;
    m_impl->mw = thermo.speciesMw(index);

    // Fill in default values
    double theta_power = std::pow(thetav, 4.0/3.0);
    size_t offset = (thermo.hasElectrons() ? 1 : 0);

    for (size_t i = 0; i < thermo.nHeavy(); ++i)
    {
        double mwi = thermo.speciesMw(i+offset);
        double mu = 1000.0*m_impl->mw*mwi/(m_impl->mw + mwi); // reduced mass in g/mol
        m_impl->a[i] = 1.16E-3*std::sqrt(mu)*theta_power;
        m_impl->b[i] = 0.015*std::pow(mu, 0.25);
    }
}

MillikanWhiteModelData::MillikanWhiteModelData(const MillikanWhiteModelData& rhs)
    : m_impl(nullptr)
{
    if (rhs.m_impl)
        m_impl = std::make_unique<Impl>(*rhs.m_impl);
}

MillikanWhiteModelData& MillikanWhiteModelData::operator=(MillikanWhiteModelData rhs)
{
    if (!rhs.m_impl)
        m_impl.reset();
    else if (!m_impl)
        m_impl = std::make_unique<Impl>(*rhs.m_impl);
    else
        *m_impl = *rhs.m_impl;
    
    return *this;
}

MillikanWhiteModelData::MillikanWhiteModelData(MillikanWhiteModelData&& rhs) noexcept = default;

MillikanWhiteModelData& MillikanWhiteModelData::operator=(MillikanWhiteModelData&& rhs) noexcept = default;

MillikanWhiteModelData::~MillikanWhiteModelData() = default;

size_t MillikanWhiteModelData::speciesIndex() const { return m_impl->index; }

double MillikanWhiteModelData::molecularWeight() const { return m_impl->mw; }

Eigen::ArrayXd& MillikanWhiteModelData::a() { return m_impl->a; }

const Eigen::ArrayXd& MillikanWhiteModelData::a() const { return m_impl->a; }

Eigen::ArrayXd& MillikanWhiteModelData::b() { return m_impl->b; }

const Eigen::ArrayXd& MillikanWhiteModelData::b() const { return m_impl->b; }

MillikanWhiteModelData& MillikanWhiteModelData::setReferenceCrossSection(double omegav)
{
    assert(omegav >= 0.0);
    m_impl->omegav = omegav;
    return *this;
}

double MillikanWhiteModelData::referenceCrossSection() const { return m_impl->omegav; }

double MillikanWhiteModelData::limitingCrossSection(const double& T) const
{
    // Park high temperature correction
    return m_impl->omegav * (2.5E9/(T*T)); // 2.5E9 = 50,000^2
}

MillikanWhiteModel::MillikanWhiteModel(const MillikanWhiteModelData& data) :
    m_data(data) 
{ }


double MillikanWhiteModel::relaxationTime(
    const Mutation::Thermodynamics::Thermodynamics& thermo) const
{
    // Millikan-White model for average relaxation time
    const Eigen::Map<const Eigen::ArrayXd> Xh(
        thermo.X()+(thermo.hasElectrons() ? 1 : 0), thermo.nHeavy());
    const double T_fac = std::pow(thermo.T(), -1.0/3.0);
    const double p_atm = thermo.P() / ONEATM;
    const double tau_mw = 
        (Xh*(m_data.a()*(T_fac - m_data.b()) - 18.42).exp()).sum() / 
        (Xh.sum()*p_atm);

    // Park correction
    const double ni = thermo.numberDensity() * thermo.X()[m_data.speciesIndex()];
    const double ci = std::sqrt(8*RU*thermo.T()/(PI*m_data.molecularWeight()));
    const double tau_park = 1.0/(ni * ci * m_data.limitingCrossSection(thermo.T()));

    return tau_mw + tau_park;
}


struct MillikanWhiteModelDB::Data
{
    const class Thermodynamics& thermo;
    XmlDocument database;
    XmlElement::const_iterator millikan_white_xml;

    Data(const class Thermodynamics& thermo, std::string file_path) :
        thermo(thermo), database(file_path)
    {
        auto& root = database.root();
        millikan_white_xml = root.findTag("Millikan-White");

        if (millikan_white_xml == root.end())
            root.parseError("Could not find Millikan-White element.");
    }
};


MillikanWhiteModelDB::MillikanWhiteModelDB(
    const class Thermodynamics& thermo
) :
    MillikanWhiteModelDB(thermo, databaseFileName("VT.xml", "transfer"))
{ }


MillikanWhiteModelDB::MillikanWhiteModelDB(
    const class Thermodynamics& thermo, std::string file_path
) :
    m_data(std::make_shared<Data>(thermo, file_path))
{ }


MillikanWhiteModel MillikanWhiteModelDB::create(
    std::string species_name, double theta)
{
    // Create default data
    MillikanWhiteModelData data(
        m_data->thermo, m_data->thermo.speciesIndex(species_name), theta);

    // Find species in database
    auto species_iter = m_data->millikan_white_xml->findTagWithAttribute(
        "vibrator", "species", species_name);
    
    // If species not in database, return default model
    if (species_iter == m_data->millikan_white_xml->end())
        return { data };
    
    // Use reference cross section if provided
    if (species_iter->hasAttribute("omegav"))
        data.setReferenceCrossSection(species_iter->getAttribute<double>("omegav"));
    
    // Loop over each heavy species in thermo
    size_t offset = (m_data->thermo.hasElectrons() ? 1 : 0);

    for (size_t i = 0; i < m_data->thermo.nHeavy(); ++i)
    {
        auto partner_iter = species_iter->findTagWithAttribute(
            "partner", "species", m_data->thermo.speciesName(i+offset));
        
        if (partner_iter != species_iter->end())
        {
            partner_iter->getAttribute("a", data.a()[i], "must provide constant a!");
            partner_iter->getAttribute("b", data.b()[i], "must provide constant b!");
        }
    }

    return { data };
}

    } // namespace Transfer
} // namespace Mutation
