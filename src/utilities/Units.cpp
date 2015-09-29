/**
 * @file Units.cpp
 *
 * @brief Implements the Units class.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdlib>

namespace Mutation {
    namespace Utilities {

/**
 * Defines all of the defined units that can be used by the user.  Units are
 * defined by first defining the 7 base SI units and then deriving all other 
 * units from these in a straight forward way.
 */
std::map<std::string, Units> _initializeUnits()
{
    std::map<std::string, Units> units;
    
    // length
    
    units["m"]    = Units(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    units["mm"]   = units["m"] / 1000.0;
    units["cm"]   = units["m"] / 100.0;
    units["km"]   = units["m"] * 1000.0;
    units["A"]    = units["m"] / 1.0e10;
    
    // mass
    units["kg"]   = Units(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    units["g"]    = units["kg"] / 1000.0;
    
    // time
    units["s"]    = Units(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    units["min"]  = units["s"] * 60.0;
    
    // thermodynamic temperature
    units["K"]    = Units(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    
    // quantity of substance
    units["mol"]  = Units(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    units["kmol"] = units["mol"] * 1000.0;
    units["molecule"] = units["mol"] * NA;
    
    // force
    units["N"]    = units["kg"] * units["m"] / units["s"] / units["s"];
    
    // energy
    units["J"]    = units["N"] * units["m"];
    units["kJ"]   = units["J"] * 1000.0;
    units["cal"]  = units["J"] * 4.184;
    units["kcal"] = units["cal"] * 1000.0;
    
    // pressure
    units["Pa"]   = units["N"] / units["m"] / units["m"];
    units["atm"]  = units["Pa"] * 101325.0;
    units["bar"]  = units["Pa"] * 100000.0;
    
    return units;
}

std::map<std::string, Units> _defined_units = _initializeUnits();

std::vector<Units> Units::split(const std::string str) {
    std::vector<std::string> tokens;
    String::tokenize(str, tokens, ",");
    std::vector<Units> units;
    for (int i = 0; i < tokens.size(); ++i)
        units.push_back(Units(tokens[i]));
    return units;
}

Units operator*(const Units& left, const double right) {
    Units units(left);    
    units.m_factor *= right;
    return units;
}

Units operator/(const Units& left, const double right) {
    Units units(left);
    units.m_factor /= right;
    return units;
}

Units operator*(const Units& left, const Units& right) {
    Units units(left);
    units.m_factor *= right.m_factor;
    for (int i = 0; i < 7; ++i)
        units.m_exponent[i] += right.m_exponent[i];
    return units;
}

Units operator/(const Units& left, const Units& right) {
    Units units(left);
    units.m_factor /= right.m_factor;
    for (int i = 0; i < 7; ++i)
        units.m_exponent[i] -= right.m_exponent[i];    
    return units;
}

Units operator^(const Units& left, const double power) {
    Units units(left);
    units.m_factor = std::pow(left.m_factor, power);
    for (int i = 0; i < 7; ++i)
        units.m_exponent[i] = 
            (left.m_exponent[i] != 0.0 ? 
                std::pow(left.m_exponent[i], power) : 0.0);
    return units;
}

Units::Units()
{
    m_factor = 1.0;
    for (int i = 0; i < 7; ++i) {
        m_exponent[i] = 0.0;
    }
}


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

Units& Units::operator=(const Units& right) {
    m_factor = right.m_factor;
    for (int i = 0; i < 7; ++i)
        m_exponent[i] = right.m_exponent[i];
    return *this;
}

Units::Units(const Units& units) {
    operator=(units);
}

Units::Units(const char cstr[])
{
    initializeFromString(std::string(cstr));
}

Units::Units(const std::string& str)
{
    initializeFromString(str);
}

void Units::initializeFromString(const std::string& str)
{
    operator=(Units());
    
    std::vector<std::string> tokens1;
    String::tokenize(str, tokens1, "/ ");
    
    std::vector<std::string> tokens2;
    String::tokenize(tokens1[0], tokens2, "- ");
    
    Units units;
    
    std::vector<std::string>::const_iterator iter;
    for (iter = tokens2.begin(); iter != tokens2.end(); ++iter) {
        if (_defined_units.find(*iter) != _defined_units.end())
            operator=(*this * _defined_units[*iter]);
        else {
            std::cerr << *iter << " is not a defined unit!" << std::endl;
            exit(1);
        }
    }
    
    if (tokens1.size() == 2) {
        std::vector<std::string> tokens3;
        String::tokenize(tokens1[1], tokens3, "- ");
        
        for (iter = tokens3.begin(); iter != tokens3.end(); ++iter) {
            if (_defined_units.find(*iter) != _defined_units.end())
                operator=(*this / _defined_units[*iter]);
            else {
                std::cerr << *iter << " is not a defined unit!" << std::endl;
                exit(1);
            }
        }
    }
}

double Units::convertToBase(const double number) {
    return number * m_factor;
}

double Units::convertTo(const double number, const Units& units) {
    for (int i = 0; i < 7; ++i)
        assert( m_exponent[i] == units.m_exponent[i] );
    return number * m_factor / units.m_factor;
}

double Units::convertTo(const double number, const std::string& units) {
    return convertTo(number, Units(units));
}

    } // namespace Utilities
} // namespace Mutation


