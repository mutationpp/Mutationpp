#ifndef DATA_GSI_RATE_MANAGER_H
#define DATA_GSI_RATE_MANAGER_H

#include "Thermodynamics.h"

#include "GSIReaction.h"
#include "SurfaceProperties.h"
#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

/**
 * Structure which stores the necessary inputs for the GSIRateManager class.
 */
struct DataGSIRateManager { 
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const SurfaceProperties& s_surf_props;
    const WallState& s_wall_state;
    const std::vector<GSIReaction*>& s_reactions;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_GSI_RATE_MANAGER_H
