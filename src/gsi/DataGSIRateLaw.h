#ifndef DATAGSIRATELAW_H
#define DATAGSIRATELAW_H 

#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

struct DataGSIRateLaw {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Utilities::IO::XmlElement& s_node_rate_law;
    const SurfaceProperties& s_surf_props;
    const std::vector<int>& s_reactants;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATAGSIRATELAW_H
