/**
 * @file ThermalConductivityAlgorithm.h
 *
 * @brief Provides ThermalConductivityAlgorithm class.
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

#ifndef TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
#define TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H

#include "CollisionDB.h"

namespace Mutation {
    namespace Transport {

/**
 * Abstract base class for all thermal conductivity algorithms which allows for 
 * self registration of concrete types.
 */
class ThermalConductivityAlgorithm
{
public:

    /// All ThermalConductivityAlgorithms must provide a constructor taking
    /// these arguements.
    typedef CollisionDB& ARGS;

    /// Returns name of this type.
    static std::string typeName() { return "ThermalConductivityAlgorithm"; }

    /**
     * Constructor.
     */
    ThermalConductivityAlgorithm(ARGS collisions)
        : m_collisions(collisions)
    { }
    
    /**
     * Destructor.
     */
    virtual ~ThermalConductivityAlgorithm() {}
    
    /**
     * Returns the mixture thermal conductivity in W/m-K using the algorithm 
     * defined in the derived type.
     */
    virtual double thermalConductivity() = 0;
    
    
    /**
     * Returns the heavy particle thermal diffusion ratios if implemented by the
     * concrete ThermalConductivityAlgorithm type.
     */
    virtual void thermalDiffusionRatios(double* const p_k)
    {
        for (int i = 0; i < m_collisions.nSpecies(); ++i)
            p_k[i] = 0.0;
    }
    

protected:
    
    CollisionDB& m_collisions;

}; // class ThermalConductivityAlgorithm

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
