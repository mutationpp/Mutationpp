#ifndef WALL_PRODUCTION_TERMS_H
#define WALL_PRODUCTION_TERMS_H

#include<eigen3/Eigen/Dense>

#include "GSIReaction.h"
#include "GSIRateManager.h"
#include "DataWallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionTerms
{
public:
    typedef const DataWallProductionTerms& ARGS;

    WallProductionTerms(ARGS args){}

	/// Returns name of this type.
	static std::string typeName() { return "WallProductionTerms"; }


    virtual ~WallProductionTerms(){}

    virtual void productionRate(Eigen::VectorXd& lv_mass_prod_rate) = 0;

    virtual const std::string& getWallProductionTermTag() const = 0;

}; // class WallProductionTerms

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // WALL_PRODUCTION_TERMS_H
