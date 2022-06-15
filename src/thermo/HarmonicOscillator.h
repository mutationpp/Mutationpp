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

#ifndef THERMO_HARMONIC_OSCILLATOR_H
#define THERMO_HARMONIC_OSCILLATOR_H

#include <memory>
#include <string>
#include <vector>

namespace Mutation {
    namespace Thermodynamics {

/**
 * Computes vibrator energy based on harmonic-oscillator model.
 */
class HarmonicOscillator
{
public:
    /// Default oscillator with zero energy.
    HarmonicOscillator() { }

    /// Construct with list of characteristic temperatures in K.
    template <typename Container>
    HarmonicOscillator(const Container&);

    /// Returns energy of vibrator in K at given temperature in K.
    double energy(double T) const;

    /// Convience function to access characteristic temperatures in K.
    const std::vector<double>& characteristicTemperatures() const {
        return m_characteristic_temps;
    }

private:
    std::vector<double> m_characteristic_temps;
};

template <typename Container>
HarmonicOscillator::HarmonicOscillator(const Container& temps) :
    m_characteristic_temps(std::begin(temps), std::end(temps))
{ }


class HarmonicOscillatorDB
{
public:
    /// Constructs the database using default file path.
    HarmonicOscillatorDB();

    /// Constructs the database using given file path.
    HarmonicOscillatorDB(std::string file_path);

    /// Creates a new harmonic oscillator associated with given species name.
    HarmonicOscillator create(std::string species_name) const;

private:
    struct Data;
    std::shared_ptr<Data> m_data;
};

    } // namespace Thermodynamics
} // namespace Mutation


#endif // THERMO_HARMONIC_OSCILLATOR_H