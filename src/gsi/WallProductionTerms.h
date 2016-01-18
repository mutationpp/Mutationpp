#ifndef WALLPRODUCTIONTERMS_H
#define WALLPRODUCTIONTERMS_H

#include<Eigen/Dense>

#include "GSIReaction.h"
#include "GSIRateManager.h"
#include "DataWallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionTerms {

public:
    typedef const DataWallProductionTerms& ARGS;

    WallProductionTerms( ARGS l_data_wall_production_terms ){ } // Check if ARGS is indeed needed
    virtual ~WallProductionTerms(){ }

    virtual void productionRate( Eigen::VectorXd& lv_mass_prod_rate ) = 0; // FIX ME!

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // WALLPRODUCTIONTERMS_H
