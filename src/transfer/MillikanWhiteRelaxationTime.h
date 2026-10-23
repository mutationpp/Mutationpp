/**
 * @file MillikanWhiteRelaxationTime.h
 *
 * @brief Relaxation time for VT source term
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


#ifndef MILLIKAN_WHITE_RELAXATION_TIME_H
#define MILLIKAN_WHITE_RELAXATION_TIME_H

#include <string>

namespace Mutation {

namespace Thermodynamics {
    class Thermodynamics;
}

namespace Transfer {

class MillikanWhiteModelData;

/**
 * @brief Interface for Millikan-White vibrational relaxation time models.
 *
 * Concrete relaxation-time models (e.g. the original Millikan-White mixing
 * rule) derive from this class and are registered.
 */
class MillikanWhiteRelaxationTime
{
public:
    using ARGS = const MillikanWhiteModelData&;

    static std::string typeName()
    {
        return "MillikanWhiteRelaxationTime";
    }

    virtual ~MillikanWhiteRelaxationTime() = default;

    /// Computes the relaxation time for the given species/mixture state.
    virtual double relaxationTime(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        const MillikanWhiteModelData& data) const = 0;

};

} // namespace Transfer
} // namespace Mutation

#endif // MILLIKAN_WHITE_RELAXATION_TIME_H
