/**
 * @file Units.h
 *
 * @brief Declares the Units class.
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


#ifndef UTILS_UNITS_H
#define UTILS_UNITS_H

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include <Eigen/Dense>

namespace Mutation {
    namespace Utilities {

class Units
{
public:

    /// Returns the map of currently defined units.
    static std::map<std::string, Units>& definedUnits();

    /// Utility function for creating a vector of Units from a comma-separated
    /// list.
    static std::vector<Units> split(const std::string str);

public:

    /// Empty constructor (dimensionless units).
    Units();

    /// Constructor taking 7 base unit dimensions.
    Units(double e1, double e2, double e3, double e4, double e5, double e6, double e7); 
    
    /// Copy constructor.
    Units(const Units& units) {
        operator=(units);
    }

    /// Constructs by parsing C-string
    Units(const char* const c_str) {
        initializeFromString(c_str);
    }

    /// Constructs by parsing string.
    Units(const std::string& str) {
        initializeFromString(str);
    }

    /**
     * Returns the factor associated with this units type.  The factor
     * corresponds to ratio between the base units and this units.  For example,
     * converting cm to m requires dividing by 100, thus the factor is 0.01.
     */
    double factor() const { return m_factor; }

    /**
     * Returns the exponents associated with the base units.  This is a 7
     * component vector corresponding to the exponents for length, mass, time,
     * electric current, thermodynamic temperature, amount of substance, and
     * luminous intensity.
     */
    Eigen::Matrix<double,7,1> exponents() const {
        Eigen::Matrix<double,7,1> exp;
        for (int i = 0; i < 7; ++i)
            exp[i] = m_exponent[i];
        return exp;
    }

    /// Converts given number from these units to appropriate base units.
    double convertToBase(const double number) {
        return m_factor * number;
    }

    /// Converts given number from these units to the given units.
    double convertTo(const double number, const Units& units);
    
    /// Returns a string representation of this Units object.
    std::string toString() const;

    /// Assignment operator.
    Units& operator=(const Units& right);

    friend Units operator*(const Units& left, const double right);
    friend Units operator/(const Units& left, const double right);
    friend Units operator*(const Units& left, const Units& right);
    friend Units operator/(const Units& left, const Units& right);
    friend Units operator^(const Units& left, const double power);

    /// Output stream operator.
    friend std::ostream& operator<< (std::ostream& out, const Units& units) {
        out << units.toString();
        return out;
    }

private:

    /// Defines the initial set of units which are defined.
    static std::map<std::string, Units> initializeUnits();

    void initializeFromString(const std::string& str);

    /// Parses a string representing a combination of units, all with POSITIVE
    /// dimension (no numerator).
    Units parseUnits(const std::string& str, const std::string& full) const;

private:

    double m_factor;
    double m_exponent[7];
};

    } // namespace Utilities
} // namespace Mutation

#endif // UTILS_UNITS_H
