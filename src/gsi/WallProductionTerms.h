#ifndef WALLPRODUCTIONTERMS_H
#define WALLPRODUCTIONTERMS_H

#include "GSIReaction.h"
#include "GSIRateManager.h"
#include "DataWallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionTerms {

public:
    typedef const DataWallProductionTerms& ARGS;

    WallProductionTerms( ARGS l_data_wall_production_terms ){ }
    virtual ~WallProductionTerms(){ }

    virtual void productionRate( Mutation::Numerics::RealVector& lv_mass_prod_rate ) = 0;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // WALLPRODUCTIONTERMS_H
