#ifndef SURFACE_BALANCE_SOLVER_H
#define SURFACE_BALANCE_SOLVER_H

#include <Eigen/Dense>

#include "NewtonSolver.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "DataSurfaceBalanceSolver.h"
#include "DiffusionVelocityCalculator.h"
#include "MassBlowingRate.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * This is the abstract class for the solution of the mass and energy balances
 * at the wall.
 */
class SurfaceBalanceSolver
{
public:
    typedef const DataSurfaceBalanceSolver& ARGS;

    /**
     * Destructor
     */
    virtual ~SurfaceBalanceSolver(){ }

    /**
     * Returns
     */
    virtual Eigen::VectorXd& computeGSIProductionRate() = 0; //  @todo const correctness // reference?

    /**
     * Purely virtual function. Temporary solution for setting up a diffusion
     * model for the surface balances.
     */
    virtual void setDiffusionModel(const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx) = 0;


    /**
     * Purely virtual function to be called in order to solve the surface balance.
     */
    virtual void solveSurfaceBalance(const Eigen::VectorXd& lv_rhoi, const Eigen::VectorXd& lv_T){}; //@BD Remove the function above and keep the one with no arguments. The state should not be passed like that, but should be retrieved by the wall state.
    virtual void solveSurfaceBalance(){ }; //@todo FIX IMPORTANT Which do I keep? Make it purely virtual

    virtual void getBprimeCondensedSpecies(std::vector<std::string>& CondensedSpecies){};
    virtual void getBprimeParameters(double & Bprime_char, std::vector<double>& CondensedMoleFrac){};

}; // class SurfaceBalanceSolver

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_BALANCE_SOLVER_H
