/**
 * @file GasFourierHeatFluxCalculator.h
 *
 * @brief Declaration of GasFourierHeatFluxCalculator class.
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


#ifndef GAS_FOURIER_HEAT_CALCULATOR_H
#define GAS_FOURIER_HEAT_CALCULATOR_H

#include <Eigen/Dense>

namespace Mutation { namespace Thermodynamics {class Thermodynamics; }}
namespace Mutation { namespace Transport {class Transport; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

/*
 * Class responsible for computing the gas to surface heat flux based on
 * Fourier's law.
 */
class GasFourierHeatFluxCalculator
{
public:
	/*
	 * Constructor
	 */
    GasFourierHeatFluxCalculator(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        Mutation::Transport::Transport& transport);

//==============================================================================
    /*
     * Default Destructor
     */
    ~GasFourierHeatFluxCalculator();

//==============================================================================
    /*
     * Function used to set up the temperature at a distance from the surface.
     *
     * @para Temperature at a distance
     * @para dx distance in m
     */
    void setGasFourierHeatFluxModel(
        const Eigen::VectorXd& v_T_edge, const double& dx);

//==============================================================================
    /*
     * Function used to compute the total fourier heat flux in
     * the gas with given the surface temperature(s) and the ones imposed at
     * a given distance from the surface. The temperature gradients
     * are computed based on a first order differentiation.
     *
     * @param Surface temperatures
     *
     */
    double computeGasFourierHeatFlux(
        const Eigen::VectorXd& v_T);

private:
    Mutation::Transport::Transport& m_transport;

    Eigen::VectorXd mv_T_edge;
    Eigen::VectorXd mv_dTdx;

    Eigen::VectorXd mv_lambda;

    double m_dx;
    bool m_is_cond_set;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GAS_FOURIER_HEAT_FLUX_CALCULATOR_H
