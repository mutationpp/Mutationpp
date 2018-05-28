#ifndef DATA_GSI_REACTION_H
#define DATA_GSI_REACTION_H

#include <string>
#include <vector>

#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "GSIRateLaw.h"
#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

/**
 * Structure which stores the necessary inputs for the GSIReaction class.
 */
struct DataGSIReaction {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Transport::Transport& s_transport;
    const SurfaceProperties& s_surf_props;
    Mutation::Utilities::IO::XmlElement::const_iterator s_iter_reaction;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_GSI_REACTION_H
