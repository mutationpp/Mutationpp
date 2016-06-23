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
    virtual Eigen::VectorXd& computeGSIProductionRates() = 0; //  @todo const correctness // reference?

    /**
     * Purely virtual function. Temporary solution for setting up a diffusion
     * model for the surface balances.
     */
    virtual void setDiffusionModel(const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx) = 0;


    /**
     * Purely virtual function to be called in order to solve the surface balance.
     */
    virtual void solveSurfaceBalance() = 0;

    virtual void getBprimeCondensedSpecies(std::vector<std::string>& CondensedSpecies){};
    virtual void getBprimeParameters(double & Bprime_char, std::vector<double>& CondensedMoleFrac){};

}; // class SurfaceBalanceSolver

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_BALANCE_SOLVER_H
