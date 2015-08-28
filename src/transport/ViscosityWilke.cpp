/**
 * @file ViscosityWilke.cpp
 *
 * @brief Provides ViscosityWilke class.
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

    ViscosityWilke(CollisionDB& collisions) 
        : ViscosityAlgorithm(collisions)
    { }

    /**
     * Returns the viscosity of the mixture in Pa-s.
     */
    double viscosity(
        double Th, double Te, double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();
        
        return wilke(
            m_collisions.etai(Th, Te, nd, p_x).tail(nh),
            m_collisions.mass().tail(nh),
            Eigen::Map<const Eigen::ArrayXd>(p_x+(ns-nh), nh));
    }
       
};

Config::ObjectProvider<ViscosityWilke, ViscosityAlgorithm> viscosityWilke("Wilke");

    } // namespace Transport
} // namespace Mutation

