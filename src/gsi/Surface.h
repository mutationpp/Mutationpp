/**
 * @file Surface.h
 *
 * @brief Declaration of the abstract Surface class.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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


#ifndef SURFACE_H
#define SURFACE_H

#include <eigen3/Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProperties;
class SurfaceState;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the Surface
 * class.
 */
struct DataSurface {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism;
    const Mutation::Utilities::IO::XmlElement& xml_feats;
    const Mutation::Utilities::IO::XmlElement& xml_surf_chem;
    const Mutation::Utilities::IO::XmlElement& xml_surf_rad;
    SurfaceState& s_surf_state;
};

//==============================================================================

/**
 * This is the abstract class for the solution of the mass and energy balances
 * at the wall.
 */
class Surface
{
public:
    typedef const DataSurface& ARGS;

//==============================================================================

	/// Returns name of this type.
	static std::string typeName() { return "Surface"; }

//==============================================================================

    /**
     * Destructor
     */
    virtual ~Surface(){ }

//==============================================================================

    /**
     * Returns
     */
    virtual void computeSurfaceReactionRates(
        Eigen::VectorXd& v_surf_chem_rates) = 0;

//==============================================================================

    /**
     * Returns the reaction rates for all the surface reactions occuring or
     * an error if no surface reactions are considered.
     */
    virtual Eigen::VectorXd computeSurfaceReactionRatesPerReaction(){
        throw LogicError()
        << "computeSurfaceReactionRatesPerReaction can be called only "
        << "when surface reactions can occur!";
    }

//==============================================================================

    /**
     * Returns the reaction rates for all the surface reactions occuring or
     * an error if no surface reactions are considered.
     */
    virtual int nSurfaceReactions(){
        throw LogicError()
        << "nSurfaceReactions can be called only "
        << "when surface reactions can occur!";
    }
//==============================================================================

    /**
     * Function for setting up the diffu
     */
    virtual void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx)
        {
            throw LogicError()
            << "setDiffusionModel can be called only when computing "
            << "the the diffusion flux to the surface are needed!";
        }

//==============================================================================

    virtual void setGasFourierHeatFluxModel(
        const Eigen::VectorXd& p_T, const double& dx)
    {
        throw LogicError()
        << "setConductiveHeatFluxModel can be called only when solving "
        << "the surface energy balance!";
    }

//==============================================================================

    virtual void setGasRadHeatFlux(const double& gas_rad_heat_flux)
    {
        throw LogicError()
        << "setGasRadHeatFlux can be called only when solving "
        << "the surface energy balance!";
    }

//==============================================================================

    /**
     * Virtual function to be called in order to solve the
     * surface balance.
     */
    virtual void solveSurfaceBalance()
    {
        throw LogicError()
        << "solveSurfaceBalance can be called only when solving "
        << "the surface energy balance!";
    }

//==============================================================================

    /**
     * Purely virtual function to be called to change the number of
     * iterations performed in the surface balance.
     */
    virtual void setIterationsSurfaceBalance(const int& iter)
    {
        throw LogicError()
        << "setIterationsSurfaceBalance can be called only when solving "
        << "the surface energy balance!";
    }

//==============================================================================

    /**
     * Purely virtual function returning the total mass blowing flux
     * due to surface and bulk phase processes.
     */
    virtual double massBlowingRate() = 0;

//==============================================================================

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_H
