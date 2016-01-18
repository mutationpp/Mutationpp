#include "DataWallProductionTerms.h"
#include "SurfaceBalanceSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

class SurfaceBalanceSolverEquilibrium : public SurfaceBalanceSolver {
public:
    SurfaceBalanceSolverEquilibrium( ARGS  l_data_surface_balance_solver ){}
    ~SurfaceBalanceSolverEquilibrium( ){}

    void computeGSIProductionRate( Eigen::VectorXd& lv_mass_prod_rate ){}

    void setDiffusionModel( const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx ){}
    void solveSurfaceBalance( const Eigen::VectorXd& lv_rhoi, const Eigen::VectorXd& lv_T ){}

};

Mutation::Utilities::Config::ObjectProvider<SurfaceBalanceSolverEquilibrium, SurfaceBalanceSolver> surface_balance_solver_equilibrium("equilibrium");
 
    } // namespace GasSurfaceInteraction
} // namespace Mutation
