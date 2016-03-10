#ifndef DATAMASSBLOWINGRATE_H
#define DATAMASSBLOWINGRATE_H

#include "Thermodynamics.h"

#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

struct DataMassBlowingRate {
	const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    WallProductionTerms& s_wall_productions_terms;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATAMASSBLOWINGRATE_H
