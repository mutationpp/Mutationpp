#ifndef DATASURFACEPROPERTIES_H
#define DATASURFACEPROPERTIES_H 

#include "Thermodynamics.h"
#include "Utilities.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

struct DataSurfaceProperties {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Utilities::IO::XmlElement& s_node_surf_props;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATASURFACEPROPERTIES_H
