#ifndef DATA_GSI_RATE_LAW_H
#define DATA_GSI_RATE_LAW_H

#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

/**
 * Structure which stores the necessary inputs for the GSIRateLaw class.
 */
struct DataGSIRateLaw {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Transport::Transport& s_transport;
    const Mutation::Utilities::IO::XmlElement& s_node_rate_law;
    const SurfaceProperties& s_surf_props;
    const std::vector<int>& s_reactants;
    const std::vector<int>& s_products;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_GSI_RATE_LAW_H
