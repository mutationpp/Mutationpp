#ifndef DATAGSIREACTION_H
#define DATAGSIREACTION_H

#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

struct DataGSIReaction {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const SurfaceProperties& s_surf_props;
    Mutation::Utilities::IO::XmlElement::const_iterator s_iter_reaction;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATAGSIREACTION_H
