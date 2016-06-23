#ifndef GAS_SURFACE_INTERACTION_H
#define GAS_SURFACE_INTERACTION_H

#include <Eigen/Dense>
#include <string>

#include "Thermodynamics.h"
#include "Transport.h"

#include "SurfaceProperties.h"
#include "WallState.h"
#include "SurfaceBalanceSolver.h" 

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 *  Interface for the computation of the Gas-Surface Interaction
 *  chemical source terms. After setting the wall state the surface
 *  chemical production rates can be returned. The solution of the
 *  surface mass balance can also be provided.
 *
 *  Currently, gamma and elementary reactions models for catalysis
 *  have been implemented and some basic models for ablation are
 *  available.
 *
 *  In the future, the solution of the energy balance at the wall
 *  will be provided and the composition of an equilibrium catalytic
 *  and ablative wall.
 *
 *  To be extended to support multiple different chemically reacting
 *  walls for a single simulation.
 *
 *  @todo Solution of the energy balance.
 *  @todo Equilibrium wall
 *  @todo Multiple chemically reacting walls.
 *
 */
class GasSurfaceInteraction
{
public:

	/**
	 * Constructor of the Gas Surface Interaction interface. Takes as inputs a
	 * reference to the thermodynamics and transport classes and the name of
	 * the gas surface interaction file. Firstly, it searches in the current
	 * directory for the file and then in /mutation++/data/gsi/
	 */
    GasSurfaceInteraction(
    	Mutation::Thermodynamics::Thermodynamics& l_thermo,
        Mutation::Transport::Transport& l_transport,
        std::string l_gsi_mechanism_file);

    /**
     * Destructor
     */
    ~GasSurfaceInteraction();

    /**
     * Function which sets the thermodynamic state of the wall. Takes as inputs
     * a pointer to the partial densities and pointer to the temperature of the
     * wall. Currently the state_variable input does not play a role.
     *
     * @param l_mass - The "mass" vector.
     * @param l_energy - The "energy" vector.
     * @param state_variable - Index representing which variable set is given in
     * the mass and energy vector.
     *
     * @todo Set the wall state according to the gas phase state model, using
     * the state_variable.
     */
    void setWallState(
        const double* const l_mass, const double* const l_energy,
		const int state_variable);

    /**
     * Function which returns the thermodynamic state of the wall. It returns
     * the values in a pointer for the partial densities and a pointer for the
     * wall temperatures. Currently, the state_variable does not play a role.
     *
     * @param l_mass - The "mass" vector.
     * @param l_energy - The "energy" vector.
     * @param state_variable - Index representing which variable set is given in
     * the mass and energy vector.
     *
     * @todo Get the wall state in the variable set according to the state
     * variable.
     */
    void getWallState(
        double* const l_mass, double* const l_energy,
        const int state_variable);

    /**
     * Return the chemical production terms of all the species in the mixture
     * due to gas surface interaction phenomena according to the current wall
     * state.
     *
     * return surface production rates in kg/m^2-s.
     */
    void surfaceProductionRates(double* const lp_wall_prod_rates);
    
    /**
     * Temporary function which set ups the diffusion model in order to compute
     * the gradient of mole fractions. Requires as input a mole fraction
     * pointer for the chemical state of the gas near the wall and a distance
     * in meters for the distance of the distance of this composition from the
     * wall.
     *
     * @param lp_mole_frac_edge     variable containing the mole fractions far
     *  from the wall used to compute the gradient for diffusion.
     * @param l_dx                  distance of the above mole fractions from
     * the wall
     */
    void setDiffusionModel(
        const double* const lp_mole_frac_edge, const double& l_dx);

    /**
     * Function to be called in order to solve the mass and energy balances at
     * the wall according to the input model. The output state is stored in the
     * wall state and can be accessed by the getWallState function.
     */
    void solveSurfaceBalance();


    /**
     * Bprime char species @BD make this function RETURN a vector of strings!
     * Chech const correctness
     */
    void getBprimeCharSpecies(std::vector<std::string>& l_species_char_names);

    /**
     * Same as above
     */
    void getBprimeSolution(
        double& l_bprime_char, std::vector<double>& lv_species_char_mass_frac);

private:
    /**
     * Helper function which locates the input file.
     */
    inline void locateGSIInputFile(std::string& l_gsi_mechanism_file);

    /**
     * Error function; wrong type of Gas Surface Interaction input file.
     */
    inline void errorWrongTypeofGSIFile(const std::string& l_gsi_root_tag);

    /**
     * Error function; invalid Gas Surface Interaction file properties.
     */
    inline void errorInvalidGSIFileProperties(const std::string& l_gsi_option);

private:

    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Transport::Transport& m_transport;

    SurfaceProperties* mp_surf_props;
    WallState* mp_wall_state;
    SurfaceBalanceSolver* mp_surf_solver;

    std::string m_gsi_mechanism;

    Eigen::VectorXd v_mass_prod_rate;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GAS_SURFACE_INTERACTION_H
