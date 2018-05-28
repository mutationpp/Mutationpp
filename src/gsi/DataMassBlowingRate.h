#ifndef DATA_MASS_BLOWING_RATE_H
#define DATA_MASS_BLOWING_RATE_H

#include "Thermodynamics.h"

#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

/**
 * Structure which stores the necessary inputs for the MassBlowingRate class.
 */
struct DataMassBlowingRate {
	const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    std::vector<WallProductionTerms*>& vs_wall_productions_terms;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_MASS_BLOWING_RATE_H
