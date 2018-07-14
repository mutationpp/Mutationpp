/**
 * @file GasSurfaceInteraction.h
 *
 * @brief Declaration of GasSurfaceInteraction class.
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


#ifndef GAS_SURFACE_INTERACTION_H
#define GAS_SURFACE_INTERACTION_H

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceBalanceSolver;
class SurfaceProperties;
class WallState;

/**
 *  Interface for the part of Mutation++ responsible for modeling Gas-Surface
 *  interaction phenomena.
 *
 *  After setting the wall state the surfac  chemical production rates can be
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
     * Function which sets the thermodynamic state of the wall. Takes as inputs
     * a pointer to the partial densities and pointer to the temperature of the
     * wall.
     *
     * @param l_mass - The "mass" vector.
     * @param l_energy - The "energy" vector.
     * @param state_variable - Index representing which variable set is given in
     * the mass and energy vector.
     */
    void setWallState(
        const double* const p_mass, const double* const p_energy,
		const int state_var);

    /**
     * Function which returns the thermodynamic state of the wall. It returns
     * the values in a pointer for the partial densities and a pointer for the
     * wall temperatures.
     *
     * @param l_mass - The "mass" vector.
     * @param l_energy - The "energy" vector.
     * @param state_variable - Index representing which variable set is given in
     * the mass and energy vector.
     */
    void getWallState(
        double* const p_mass, double* const p_energy,
        const int state_var);

    /**
     * Return the chemical production terms of all the species in the mixture
     * due to gas surface interaction phenomena according to the current wall
     * state.
     *
     * return surface production rates in kg/m^2-s.
     */
    void surfaceProductionRates(double* const p_wall_prod_rates);
    
    /**
     * Temporary function which set ups the diffusion model in order to compute
     * the gradient of mole fractions. Requires as input a mole fraction
     * pointer for the chemical state of the gas near the wall and a distance
     * in meters for the distance of the distance of this composition from the
     * wall.
     *
     * @param mole_frac_edge   mole fractions at a distance from the wall
     *                         used to compute the gradient for diffusion
     * @param dx               distance from the wall in m
     */
    void setDiffusionModel(
        const double* const p_mole_frac_edge, const double& dx);

    /*
     *  Function for the energy balance at the wall. Not used yet.
     */
    void setConductiveHeatFluxModel(
        const double* const p_T_edge, const double& dx);

    /**
     * Function to be called in order to solve the mass and energy balances at
     * the wall according to the input model. The output state is stored in the
     * wall state and can be accessed by the getWallState function.
     */
    void solveSurfaceBalance();

    /**
     * Function which return the total mass blowing flux.
     *
     * @param mdot on return mass blowing flux kg/(m^2-s)
     */
    void getMassBlowingRate(double& mdot);

    /**
     * Function for the wall in equilibrium. Not used yet.
     */
    void getBprimeCharSpecies(std::vector<std::string>& v_species_char_names);

    /**
     * Function for the wall in equilibrium. Not used yet.
     */
    void getBprimeSolution(
        double& bprime_char, std::vector<double>& v_species_char_mass_frac);

private:
    /**
     * Error function; wrong type of Gas Surface Interaction input file.
     */
    inline void errorWrongTypeofGSIFile(const std::string& gsi_root_tag);

    /**
     * Error function; invalid Gas Surface Interaction file properties.
     */
    inline void errorInvalidGSIFileProperties(const std::string& gsi_option);

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Transport::Transport& m_transport;

    SurfaceProperties* mp_surf_props;
    WallState* mp_wall_state;
    SurfaceBalanceSolver* mp_surf_solver;

    std::string m_gsi_mechanism;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GAS_SURFACE_INTERACTION_H
