/**
 * @file DiffusionVelocityCalculator.h
 *
 * @brief Declaration of DiffustionVelocityCalculator class.
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


#ifndef DIFFUSION_VELOCITY_CALCULATOR_H
#define DIFFUSION_VELOCITY_CALCULATOR_H

#include <Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

/*
 * Class responsible for computing the diffusion velocity in the gas.
 */
class DiffusionVelocityCalculator
{
public:
	/*
	 * Constructor of the class. It takes as arguments a copy of the
	 * thermodynamics class and a copy of the transport class in order
	 * to call the Stefan-Maxwell method.
	 */
    DiffusionVelocityCalculator(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        Mutation::Transport::Transport& transport);

//==============================================================================
    /*
     * Default Destructor
     */
    ~DiffusionVelocityCalculator();

//==============================================================================
    /*
     * Function used to set up the mole fractions at a distance from the
     * surface.
     *
     * @para mole_frac_edge mole fractions of the species at a distance
     * @para dx distance in m
     */
    void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx);

//==============================================================================
    /*
     * Function used to compute the diffusion velocities with the given
     * mole fraction at the surface the ones imposed at a given distance from
     * the surface after the diffusion model has been set. The mole fraction
     * gradients are computed with simple first order differentiation.
     *
     * @param mole_frac mole fractions of the species on the surface
     * @param on return diffusion velocities in m/s
     *
     */
    void computeDiffusionVelocities(
        const Eigen::VectorXd& v_mole_frac,
        Eigen::VectorXd& v_diff_velocities);

private:
    Mutation::Transport::Transport& m_transport;

    Eigen::VectorXd mv_mole_frac_edge;
    Eigen::VectorXd mv_dxidx;

    double m_dx;
    bool m_is_diff_set;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSION_VELOCITY_CALCULATOR_H
