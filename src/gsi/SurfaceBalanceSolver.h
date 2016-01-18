#ifndef SURFACEBALANCESOLVER_H
#define SURFACEBALANCESOLVER_H

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

class SurfaceBalanceSolver {
public:
    typedef const DataSurfaceBalanceSolver& ARGS;

    virtual ~SurfaceBalanceSolver(){ }

    virtual void computeGSIProductionRate( Eigen::VectorXd& lv_mass_prod_rate ) = 0;

    virtual void setDiffusionModel( const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx ) = 0;
    virtual void solveSurfaceBalance( const Eigen::VectorXd& lv_rhoi, const Eigen::VectorXd& lv_T ) = 0;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACEBALANCESOLVER_H
