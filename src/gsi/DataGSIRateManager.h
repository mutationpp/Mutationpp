#ifndef DATAGSIRATEMANAGER_H
#define DATAGSIRATEMANAGER_H

#include "Thermodynamics.h"

#include "GSIReaction.h"
#include "SurfaceDescription.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

struct DataGSIRateManager { 
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const SurfaceDescription& s_surf_descr;
    const std::vector<GSIReaction*>& s_reactions;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATAGSIRATEMANAGER_H
