/**
 * @file GasSurfaceInteraction.h
 *
 * @brief Declaration of GasSurfaceInteraction class.
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


#ifndef GAS_SURFACE_INTERACTION_H
#define GAS_SURFACE_INTERACTION_H

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class Surface;
class SurfaceState;

/**
 *  Interface for the part of Mutation++ responsible for modeling Gas-Surface
 *  interaction phenomena.
 *
 *  After setting the surface state the surfac  chemical production rates can be
 *  returned. The solution of the surface mass balance can also be provided.
 *
 *  Currently, gamma models for catalysis and ablation are available.
 */

class GasSurfaceInteraction
{
public:

	/**
	 * Constructor of the GasSurfaceInteraction class. Takes as inputs a
	 * reference to the thermodynamics and transport classes and the name of
	 * the gas surface interaction file.
	 */
    GasSurfaceInteraction(
    	Mutation::Thermodynamics::Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        std::string gsi_mechanism_file);

    /**
     * Destructor
     */
    ~GasSurfaceInteraction();

    /**
     * Function which sets the thermodynamic state of the surface. Takes as inputs
     * a pointer to the partial densities and pointer to the temperature of the
     * surface.
     *
     * @param l_mass - The "mass" vector.
     * @param l_energy - The "energy" vector.
     * @param state_variable - Index representing which variable set is given in
     * the mass and energy vector.
     */
    void setSurfaceState(
        const double* const p_mass, const double* const p_energy,
		const int state_var);

    /**
     * Function which returns the thermodynamic state of the surface. It returns
     * the values in a pointer for the partial densities and a pointer for the
     * surface temperatures.
     *
     * @param l_mass - The "mass" vector.
     * @param l_energy - The "energy" vector.
     * @param state_variable - Index representing which variable set is given in
     * the mass and energy vector.
     */
    void getSurfaceState(
        double* const p_mass, double* const p_energy,
        const int state_var);

    /**
     * Return the chemical production terms of all the species in the mixture
     * due to gas surface interaction phenomena according to the current surface
     * state.
     *
     * return surface production rates in kg/m^2-s.
     */
    void surfaceReactionRates(double* const p_surface_reac_rates);

    /**
     * Returns the chemical production rates of all reactions
     * due to gas surface interaction phenomena according to the current surface
     * state.
     *
     * return surface production rates for each reaction.
     */
    void surfaceReactionRatesPerReaction(double* const p_reaction_rates);

    /**
     * Returns the number of surface reactions occuring during
     * Gas-Surface Interaction phenomena.
     */
    int nSurfaceReactions();

    /**
     * Function which set ups the diffusion model in order to compute
     * the gradient of mole fractions. Requires as input a mole fraction
     * pointer for the chemical state of the gas near the surface and a distance
     * these mole fractions are computed in meters.
     *
     * @param p_mole_frac_edge   mole fractions at a distance from the surface
     *                           used to compute the gradient for diffusion
     * @param dx                 distance from the surface in m
     */
    void setDiffusionModel(
        const double* const p_mole_frac_edge, const double& dx);

     /**
     * Function which set ups the conductive heat flux model necessary for the
     * energy balance. Requires as input the temperature vector at a distance
     * from the surface and the distance in meters.
     *
     * @param p_T_edge         temperatures at a distance from the surface
     *                         used to compute the gradient for diffusion
     * @param dx               distance from the surface in m
     */
    void setGasFourierHeatFluxModel(
        const double* const p_T_edge, const double& dx);

    /*
     * Function for the energy balance at the surface. Works only when the
     * surface_feature gas radiation is on. Otherwise the Tenv can be imposed
     * to compute the far field incoming radiation.
     */
    void setGasRadHeatFlux(const double* const m_gas_rad_heat_flux);

    /*
     * Function for the energy balance at the surface. Works only when the
     * surface_feature solid_conduction is set to input. Otherwise conduction is
     * computed based on the steady state model or it is neglected.
     */
    void setSolidCondHeatFlux(const double* const m_solid_heat_flux){}


    /**
     * Function to be called in order to solve the mass and energy balances at
     * the surface according to the input model. The output state is stored in the
     * surface state and can be accessed by the getSurfaceState function.
     */
    void solveSurfaceBalance();

    /**
     * Function which allows to change the number of iterations when solving
     * the surface balance at the interface. The default value is 5.
     */
    void setIterationsSurfaceBalance(const int& iter);

    /**
     * Function which return the total mass blowing flux.
     *
     * @param mdot on return mass blowing flux kg/(m^2-s)
     */
    void getMassBlowingRate(double& mdot);

private:
    /**
     * Error function; wrong type of Gas Surface Interaction input file.
     */
    inline void errorWrongTypeofGSIFile(const std::string& gsi_root_tag);

    /**
     * Error function; invalid Gas Surface Interaction file properties.
     */
    inline void errorInvalidGSIFileProperties(const std::string& gsi_option);

    /**
     * Error function; No solid properties were provided.
     */
    inline void errorSolidPropertiesNotProvided(
        const std::string& error_steady_state);

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Transport::Transport& m_transport;

    SurfaceState* mp_surf_state;
    Surface* mp_surf;

    std::string m_gsi_mechanism;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GAS_SURFACE_INTERACTION_H
