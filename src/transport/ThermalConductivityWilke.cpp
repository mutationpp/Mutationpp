/**
 * @file ThermalConductivityWilke.cpp
 *
 * @brief Provides ThermalConductivityWilke class.
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

#include "ThermalConductivityAlgorithm.h"
#include "Wilke.h"
#include "Constants.h"
#include "Numerics.h"
#include "CollisionDB.h"
#include "Utilities.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

/**
 * Computes the translational thermal conductivity of a mixture using the Wilke 
 * mixture rule and Eucken's relation.
 */
class ThermalConductivityWilke 
    : public ThermalConductivityAlgorithm, public Wilke
{
public:
 
    ThermalConductivityWilke(ThermalConductivityAlgorithm::ARGS arguments)
        : ThermalConductivityAlgorithm(arguments)
    { }
    
    double thermalConductivity(
        double Th, double Te, double nd, const double* const p_x) 
    {
        const size_t ns = m_collisions.nSpecies();
        const RealVector Mw = m_collisions.mass() * NA;
        
        return wilke(
            3.75 * RU * m_collisions.etai(Th, Te, nd, p_x) / Mw, 
            m_collisions.mass(), asVector(p_x, ns));
    }
};

// Register the algorithm
Config::ObjectProvider<ThermalConductivityWilke, ThermalConductivityAlgorithm>
    lambdaWilke("Wilke");

    } // namespace Transport
} // namespace Mutation

