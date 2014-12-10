/**
 * @file NumConst.h
 *
 * @brief Implements NumConst class.
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

#ifndef NUMERICS_NUMCONST_H
#define NUMERICS_NUMCONST_H

#include <limits>

namespace Mutation {
    namespace Numerics {


/**
 * Stores numeric constants for a template type which are used throughout the
 * algorithms in the numerics directory.
 */
template <typename T>
class NumConst
{
public:
    static const T zero;
    static const T one;
    static const T two;
    static const T eps; ///< Machine round-off
};

template <typename T>
const T NumConst<T>::zero = static_cast<T>(0);

template <typename T>
const T NumConst<T>::one  = static_cast<T>(1);

template <typename T>
const T NumConst<T>::two  = static_cast<T>(2);

template <typename T>
const T NumConst<T>::eps  = std::numeric_limits<T>::epsilon();


    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_NUMCONST_H
