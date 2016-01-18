#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"

#include <Eigen/Dense>

#include "NewtonSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerFRC : public GSIRateManager, public Mutation::Numerics::NewtonSolver<Eigen::VectorXd, GSIRateManagerFRC> {

public:
    GSIRateManagerFRC( DataGSIRateManager l_data_rate_manager ) 
                     : GSIRateManager( l_data_rate_manager ){ }

//=============================================================================

    ~GSIRateManagerFRC(){ }

//=============================================================================

    void computeRate( Eigen::VectorXd& lv_mass_prod_rate ){
        computeSteadyStateWallSpecies();
        // computeGasSpeciesProductionRate( );
    }

//=============================================================================

private:
    void computeSteadyStateWallSpecies(){ }

};

Mutation::Utilities::Config::ObjectProvider<GSIRateManagerFRC, GSIRateManager> gsi_rate_manager_frc("frc");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
