/**
 * @file Units.cpp
 *
 * @brief Implements the Units class.
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

#include "Utilities.h" // includes Units.h
#include "Constants.h"
#include "Errors.h"

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <algorithm>

namespace Mutation {
    namespace Utilities {

//==============================================================================

std::map<std::string, Units>& Units::definedUnits()
{
    static std::map<std::string, Units> defined_units = initializeUnits();
    return defined_units;
}

//==============================================================================

std::map<std::string, Units> Units::initializeUnits()
{
    std::map<std::string, Units> units;
    
    // length
    units["m"]    = Units(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    units["mm"]   = units["m"] / 1000.0;
    units["cm"]   = units["m"] / 100.0;
    units["km"]   = units["m"] * 1000.0;
    units["Å"]    = units["m"] / 1.0e10; // Å is \u00C5 in unicode
    
    // mass
    units["kg"]   = Units(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    units["g"]    = units["kg"] / 1000.0;
    
    // time
    units["s"]    = Units(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    units["min"]  = units["s"] * 60.0;
    
    // electric current
    units["A"]    = Units(0, 0, 0, 1, 0, 0, 0);

    // thermodynamic temperature
    units["K"]    = Units(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    
    // quantity of substance
    units["mol"]  = Units(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    units["kmol"] = units["mol"] * 1000.0;
    units["molecule"] = units["mol"] / NA;
    
    // luminous intensity
    units["cd"] = Units(0, 0, 0, 0, 0, 0, 1);

    // force
    units["N"]    = units["kg"] * units["m"] / units["s"] / units["s"];
    
    // energy
    units["J"]    = units["N"] * units["m"];
    units["kJ"]   = units["J"] * 1000.0;
    units["cal"]  = units["J"] * 4.184;
    units["kcal"] = units["cal"] * 1000.0;
    units["eV"]   = units["J"] / 6.241509e18;
    units["meV"]  = units["eV"] / 1000.0;
    
    // pressure
    units["Pa"]   = units["N"] / units["m"] / units["m"];
    units["atm"]  = units["Pa"] * 101325.0;
    units["bar"]  = units["Pa"] * 100000.0;
    
    return units;
}

//==============================================================================

std::vector<Units> Units::split(const std::string str)
{
    std::vector<std::string> tokens;
    String::tokenize(str, tokens, ",");
    std::vector<Units> units;
    for (int i = 0; i < tokens.size(); ++i)
        units.push_back(Units(tokens[i]));
    return units;
}

//==============================================================================

Units::Units()
{
    m_factor = 1.0;
    for (int i = 0; i < 7; ++i) {
        m_exponent[i] = 0.0;
    }
}

//==============================================================================

Units::Units(
    double e1, double e2, double e3, double e4, double e5, double e6, double e7)
    : m_factor(1.0)
{ 
    m_exponent[0] = e1;
    m_exponent[1] = e2;
    m_exponent[2] = e3;
    m_exponent[3] = e4;
    m_exponent[4] = e5;
    m_exponent[5] = e6;
    m_exponent[6] = e7;
}

//==============================================================================

Units& Units::operator=(const Units& right)
{
    m_factor = right.m_factor;
    for (int i = 0; i < 7; ++i)
        m_exponent[i] = right.m_exponent[i];
    return *this;
}

//==============================================================================

void Units::initializeFromString(const std::string& to_parse)
{
    // Make sure string is to not empty or has '/' at front or back
    std::string str = String::trim(to_parse);
    if (str.size() == 0 || str[0] == '/' || str[str.size()-1] == '/') {
        throw InvalidInputError("units", to_parse)
            << "String is empty or has empty numerator or denominator.";
    }

    // Tokenize around '/'
    std::vector<std::string> tokens;
    String::tokenize(str, tokens, "/");
    
    // After above, the only possibilities are tokens.size() > 0
    if (tokens.size() == 1) {
        operator=(parseUnits(tokens[0], str));
    } else if (tokens.size() == 2) {
        operator=(parseUnits(tokens[0], str) / parseUnits(tokens[1], str));
    } else {
        throw InvalidInputError("units", str)
            << "Can only include the \"/\" character once.";
    }
}

//==============================================================================

Units Units::parseUnits(const std::string& str, const std::string& full) const
{
    // Tokenize the string
    std::vector<std::string> tokens;
    String::tokenize(str, tokens, " .-");

    Units units;

    std::vector<std::string>::const_iterator it;
    for (it = tokens.begin(); it != tokens.end(); ++it) {
        if (definedUnits().find(*it) != definedUnits().end())
            units = units * definedUnits()[*it];
        else {
            InvalidInputError error("units", full);
            error << "\"" << *it << "\" is not a defined unit. Available "
                  << "units are:";
            std::map<std::string, Units>::const_iterator it =
                definedUnits().begin();
            for ( ; it != definedUnits().end(); ++it) {
                error << "\n  " << it->first;
            }
            throw error;
        }
    }

    return units;
}

//==============================================================================

double Units::convertTo(const double number, const Units& units)
{
    // Make sure both units represent the same quantity
    for (int i = 0; i < 7; ++i) {
        if (m_exponent[i] == units.m_exponent[i])
            continue;

        throw Error("invalid unit conversion")
            << "\"" << (*this) << "\" cannot be converted to \""
            << units << "\" because they are not composed of the same base "
            << "units.";
    }

    return number * m_factor / units.m_factor;
}

//==============================================================================

std::string Units::toString() const
{
    const char* base_units [7] = {"m", "kg", "s", "A", "K", "mol", "cd"};

    std::stringstream ss;

    // Create the numerator
    bool has_denominator = false;
    for (int i = 0; i < 7; ++i) {
        if (m_exponent[i] > 0) {
            ss << base_units[i];
            if (m_exponent[i] != 1)
                ss << "^" << m_exponent[i];
            ss << " ";
        } else if (m_exponent[i] < 0)
            has_denominator = true;
    }

    // Create the denominator
    if (has_denominator) {
        ss << "/";
        for (int i = 0; i < 7; ++i) {
            if (m_exponent[i] < 0) {
                ss << " " << base_units[i];
                if (m_exponent[i] != -1)
                    ss << "^" << std::abs(m_exponent[i]);
            }
        }
    }

    return ss.str();
}

//==============================================================================

Units operator*(const Units& left, const double right) {
    Units units(left);
    units.m_factor *= right;
    return units;
}

//==============================================================================

Units operator/(const Units& left, const double right) {
    Units units(left);
    units.m_factor /= right;
    return units;
}

//==============================================================================

Units operator*(const Units& left, const Units& right) {
    Units units(left);
    units.m_factor *= right.m_factor;
    for (int i = 0; i < 7; ++i)
        units.m_exponent[i] += right.m_exponent[i];
    return units;
}

//==============================================================================

Units operator/(const Units& left, const Units& right) {
    Units units(left);
    units.m_factor /= right.m_factor;
    for (int i = 0; i < 7; ++i)
        units.m_exponent[i] -= right.m_exponent[i];
    return units;
}

//==============================================================================

Units operator^(const Units& left, const double power) {
    Units units(left);
    units.m_factor = std::pow(left.m_factor, power);
    for (int i = 0; i < 7; ++i)
        units.m_exponent[i] *= power;
    return units;
}

//==============================================================================

    } // namespace Utilities
} // namespace Mutation


