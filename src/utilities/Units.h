/**
 * @file Units.h
 *
 * @brief Declares the Units class.
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


#ifndef UTILS_UNITS_H
#define UTILS_UNITS_H

#include <vector>
#include <string>

namespace Mutation {
    namespace Utilities {

class Units
{
public:

    static std::vector<Units> split(const std::string str);

public:

    Units();
    Units(const char []);
    Units(const std::string& str);
    //Units(std::initializer_list<double> exponents);
    Units(double e1, double e2, double e3, double e4, double e5, double e6, double e7); 
    Units(const Units& units);
    
    double convertToBase(const double number);
    double convertTo(const double number, const Units& units);
    double convertTo(const double number, const std::string& units);
    
    Units& operator=(const Units& right);
    
    friend Units operator*(const Units& left, const double right);
    friend Units operator/(const Units& left, const double right);
    friend Units operator*(const Units& left, const Units& right);
    friend Units operator/(const Units& left, const Units& right);
    friend Units operator^(const Units& left, const double power);

private:

    void initializeFromString(const std::string& str);

private:

    double m_factor;
    double m_exponent[7];
};

    } // namespace Utilities
} // namespace Mutation

#endif // UTILS_UNITS_H
