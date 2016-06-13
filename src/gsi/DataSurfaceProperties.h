#ifndef DATA_SURFACE_PROPERTIES_H
#define DATA_SURFACE_PROPERTIES_H

#include "Thermodynamics.h"
#include "Utilities.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

/**
 * Structure which stores the necessary inputs for the SurfaceProperties class.
 */
struct DataSurfaceProperties {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Utilities::IO::XmlElement& s_node_surf_props;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_SURFACE_PROPERTIES_H
