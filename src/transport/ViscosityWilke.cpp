/**
 * @file ViscosityWilke.cpp
 *
 * @brief Provides ViscosityWilke class.
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

#include "AutoRegistration.h"
#include "CollisionDB.h"
#include "Thermodynamics.h"
#include "ViscosityAlgorithm.h"
#include "Wilke.h"

#include <Eigen/Dense>

using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

/**
 * Computes the viscosity of a mixture in Pa-s using the Wilke mixture rule.
 */
class ViscosityWilke : public ViscosityAlgorithm, public Wilke
{
public:

    ViscosityWilke(ViscosityAlgorithm::ARGS collisions)
        : ViscosityAlgorithm(collisions)
    { }

    /// Returns the viscosity of the mixture in Pa-s.
    double viscosity()
    {
        const int nh = m_collisions.nHeavy();
        const int k  = m_collisions.nSpecies()-nh;

        return wilke(
            m_collisions.etai(),
            m_collisions.mass().tail(nh),
            Eigen::Map<const Eigen::ArrayXd>(m_collisions.thermo().X()+k, nh));
    }

};

// Register this algorithm
Config::ObjectProvider<ViscosityWilke, ViscosityAlgorithm> visc_wilke("Wilke");

    } // namespace Transport
} // namespace Mutation

