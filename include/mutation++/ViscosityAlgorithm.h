/**
 * @file ViscosityAlgorithm.h
 *
 * @brief Provides ViscosityAlgorithm class.
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

#ifndef TRANSPORT_VISCOSITY_ALGORITHM_H
#define TRANSPORT_VISCOSITY_ALGORITHM_H

namespace Mutation {
    namespace Transport {

class CollisionDB;

/**
 * Abstract base class for all viscosity algorithms which allows for self 
 * registration of concrete types.
 */
class ViscosityAlgorithm
{
public:
    
    // Required for self registering viscosity algorithms
    typedef CollisionDB& ARGS;
    
    /// Returns name of this type.
    static std::string typeName() { return "ViscosityAlgorithm"; }

    /// Constructor taking a CollisionDB object.
    ViscosityAlgorithm(ARGS collisions)
        : m_collisions(collisions)
    { }
    
    /// Destructor.
    virtual ~ViscosityAlgorithm() { }
    
    /// Returns the mixture viscosity in Pa-s.
    virtual double viscosity() = 0;

protected:

    CollisionDB& m_collisions;
    
}; // class ViscosityAlgorithm

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_VISCOSITY_ALGORITHM_H
