/*
 * Copyright 2014-2021 von Karman Institute for Fluid Dynamics (VKI)
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

#include "HarmonicOscillator.h"
#include "Utilities.h"

#include <cassert>
#include <cmath>
#include <iterator>
#include <sstream>

namespace Mutation {
    namespace Thermodynamics {

double HarmonicOscillator::energy(double T) const
{
    assert(T > 0.0);
    
    double energy = 0.0;
    for (auto theta: m_characteristic_temps)
        energy += theta / (std::exp(theta/T) - 1.0);
    return energy;
}


struct HarmonicOscillatorDB::Data
{
    Utilities::IO::XmlDocument document;

    Data(std::string file_path) : document(file_path) { }
};

HarmonicOscillatorDB::HarmonicOscillatorDB() :
    HarmonicOscillatorDB(Utilities::databaseFileName("species.xml", "thermo")) 
{ }

HarmonicOscillatorDB::HarmonicOscillatorDB(std::string file_path) :
    m_data(std::make_shared<Data>(file_path))
{ }

HarmonicOscillator HarmonicOscillatorDB::create(
    std::string species_name) const
{
    auto& root = m_data->document.root();
    auto iter = root.findTagWithAttribute("species", "name", species_name);
    
    if (iter == root.end())
        return {};
    
    // Find thermodynamic data
    auto thermo_iter = iter->findTagWithAttribute("thermodynamics", "type", "RRHO");
    if (thermo_iter == iter->end())
        return {};

    // Find vibrational temperature data
    iter = thermo_iter->findTag("vibrational_temperatures");
    if (iter == thermo_iter->end())
        return {};

    // Split into individual temperatures
    std::istringstream ss(iter->text());
    std::vector<double> temps;

    std::copy(
        std::istream_iterator<double>(ss), std::istream_iterator<double>(),
        std::back_inserter(temps)
    );

    return { temps };
}

    } // namespace Thermodynamics
} // namespace Mutation