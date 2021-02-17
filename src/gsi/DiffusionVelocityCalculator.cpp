/**
 * @file DiffusionVelocityCalculator.cpp
 *
 * @brief Class which computes the diffusion velocities needed by
 *        the surface balances.
 */

/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
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


#include "Errors.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "DiffusionVelocityCalculator.h"

using namespace Eigen;

namespace Mutation {
    namespace GasSurfaceInteraction {

DiffusionVelocityCalculator::DiffusionVelocityCalculator(
    const Mutation::Thermodynamics::Thermodynamics& thermo,
    Mutation::Transport::Transport& transport)
        : mv_dxidx(thermo.nSpecies()),
          mv_mole_frac_edge(thermo.nSpecies()),
          m_transport(transport),
          m_dx(0.),
          m_is_diff_set(false)
{ }

//==============================================================================

DiffusionVelocityCalculator::~DiffusionVelocityCalculator(){}

//==============================================================================

void DiffusionVelocityCalculator::setDiffusionModel(
    const VectorXd& v_mole_frac_edge, const double& dx)
{
    mv_mole_frac_edge = v_mole_frac_edge;

    if (dx <= 0.) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::setDiffusionModel() with a "
        << "distance less or equal to zero. The distance dx should always be "
           "positive.";
    }
    m_dx = dx;

    m_is_diff_set = true;
}

//==============================================================================

void DiffusionVelocityCalculator::computeDiffusionVelocities(
    const VectorXd& v_mole_frac,
    VectorXd& v_diff_velocities)
{
    if (!m_is_diff_set) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::computeDrivingForces() before "
        << "calling DiffusionVelocityCalculator::setDiffusionCalculator().";
    }
    mv_dxidx = (v_mole_frac - mv_mole_frac_edge)/m_dx;

    double electric_field = 0.E0;
    m_transport.stefanMaxwell(
        mv_dxidx.data(),
        v_diff_velocities.data(),
        electric_field);
}

    } // namespace GasSurfaceInteraction
} // namespace Mutation
