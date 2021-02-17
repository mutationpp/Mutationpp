/**
 * @file Constants.h
 *
 * @brief Defines several constants which are used throughout Mutation++.
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


#ifndef GENERAL_CONSTANTS_H
#define GENERAL_CONSTANTS_H

#include <cmath>

namespace Mutation {

/// \defgroup constants
/// @{

/// \defgroup math Mathematical constants
/// @{
const double PI    = 4.0 * std::atan(1.0); //!< \f$\pi\f$
const double TWOPI = 2.0 * PI;             //!< \f$2\pi\f$
const double SQRT2  = std::sqrt(2.0);      //!< \f$\sqrt{2}\f$
/// @}

/// \defgroup physical Physical constants (all units in SI)
/// @{
const double NA     = 6.0221415E23;    //!< Avagadro's number (molecule/mol)
const double KB     = 1.3806503E-23;   //!< Boltzmann's constant (J/molecule-K)
const double RU     = NA * KB;         //!< Universal Gas constant (J/mole-K)
const double HP     = 6.626068E-34;    //!< Planck's constant (J-s)
const double MU0    = PI * 4.0E-7;     //!< Magnetic constant (H/m)
const double C0     = 299792458.0;     //!< Speed of light in vacuum (m/s)
const double EPS0   = 1.0/(MU0*C0*C0); //!< Vacuum permittivity (F/m)
const double QE     = 1.602176565E-19; //!< Elementary positive charge (C)
const double ONEATM = 101325.0;        //!< 1 atm in Pa
/// @}

/// @}

} // namespace Mutation

#endif // GENERAL_CONSTANTS_H
