/**
 * @file GasFourierHeatFluxCalculator.cpp
 *
 * @brief Class which computes the gas heat flux needed by
 *        the surface energy balances.
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

#include "GasFourierHeatFluxCalculator.h"

using namespace Eigen;

namespace Mutation {
    namespace GasSurfaceInteraction {

GasFourierHeatFluxCalculator::GasFourierHeatFluxCalculator(
    const Mutation::Thermodynamics::Thermodynamics& thermo,
    Mutation::Transport::Transport& transport)
        : mv_dTdx(thermo.nEnergyEqns()),
          mv_T_edge(thermo.nEnergyEqns()),
          mv_lambda(thermo.nEnergyEqns()),
          m_transport(transport),
          m_dx(0.),
          m_is_cond_set(false)
{ }

//==============================================================================

GasFourierHeatFluxCalculator::~GasFourierHeatFluxCalculator(){}

//==============================================================================

void GasFourierHeatFluxCalculator::setGasFourierHeatFluxModel(
    const VectorXd& v_T_edge, const double& dx)
{
    mv_T_edge = v_T_edge;

    if (dx <= 0.) {
    	throw LogicError()
        << "Calling GasFourierHeatFluxCalculator::setGasFourierHeatFluxModel()"
        << " with a distance less or equal to zero. The distance dx should "
        << "always be positive.";
    }
    m_dx = dx;

    m_is_cond_set = true;
}

//==============================================================================

double GasFourierHeatFluxCalculator::computeGasFourierHeatFlux(
    const VectorXd& v_T)
{
    if (!m_is_cond_set) {
    	throw LogicError()
        << "Calling GasFourierHeatFluxCalculator::"
        << "HeatFluxCalculator() before "
        << "calling GasFourierHeatFluxCalculator::"
        << "setGasFourierHeatFluxModel().";
    }

    mv_dTdx = (v_T - mv_T_edge)/m_dx;
    m_transport.frozenThermalConductivityVector(mv_lambda.data());
    return -mv_lambda.dot(mv_dTdx);
}

    } // namespace GasSurfaceInteraction
} // namespace Mutation