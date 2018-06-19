#ifndef SURFACE_BALANCE_SOLVER_H
#define SURFACE_BALANCE_SOLVER_H

#include <eigen3/Eigen/Dense>

namespace Mutation {
    namespace GasSurfaceInteraction {

class Thermodynamics;
class Transport;
class XmlElement;

class SurfaceProperties;
class WallState;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the SurfaceBalanceSolver
 * class.
 */
struct DataSurfaceBalanceSolver {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism;
    const Mutation::Utilities::IO::XmlElement& s_node_diff_model;
    const Mutation::Utilities::IO::XmlElement& s_node_prod_terms;
    SurfaceProperties& s_surf_props;
    WallState& s_wall_state;
};

//==============================================================================

/**
 * This is the abstract class for the solution of the mass and energy balances
 * at the wall.
 */
class SurfaceBalanceSolver
{
public:
    typedef const DataSurfaceBalanceSolver& ARGS;

	/// Returns name of this type.
	static std::string typeName() { return "SurfaceBalanceSolver"; }

    /**
     * Destructor
     */
    virtual ~SurfaceBalanceSolver(){ }

    /**
     * Returns
     */
    virtual Eigen::VectorXd computeGSIProductionRates() = 0;

    /**
     * Function for setting up the diffu
     */
    virtual void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx) = 0;

    virtual void setConductiveHeatFluxModel(
        const Eigen::VectorXd& p_T, const double& dx)
    {
        throw LogicError()
        << "setConductiveHeatFluxModel can be called only when solving "
        << "the surface energy balance!";
    }

    /**
     * Purely virtual function to be called in order to solve the
     * surface balance.
     */
    virtual void solveSurfaceBalance() = 0;

    /**
     * Purely virtual function returning the total mass blowing flux
     * due to surface and bulk phase processes.
     */
    virtual double massBlowingRate() = 0;

    virtual void getBprimeCondensedSpecies(
        std::vector<std::string>& CondensedSpecies)
    {
        throw LogicError()
        << "getBprimeCondensedSpecies can be called only for a "
        << "surface in Equilibrium!";
    }

    virtual void getBprimeParameters(
        double & Bprime_char, std::vector<double>& CondensedMoleFrac)
    {
        throw LogicError()
        << "getBprimeParameters can be called only for a "
        << "surface in Equilibrium!";
    }

}; // class SurfaceBalanceSolver

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_BALANCE_SOLVER_H
