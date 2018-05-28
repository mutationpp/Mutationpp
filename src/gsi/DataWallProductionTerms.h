#ifndef DATA_WALL_PRODUCTION_TERMS_H
#define DATA_WALL_PRODUCTION_TERMS_H

#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "SurfaceProperties.h"
#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

class WallProductionTerms;

/**
 * Structure which stores the necessary inputs for the WallProductionTerme class.
 */
struct DataWallProductionTerms {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism; 
    const Mutation::Utilities::IO::XmlElement& s_node_prod_terms;
    const SurfaceProperties& s_surf_props;
    const WallState& s_wall_state;
    std::vector<WallProductionTerms*>* sp_surf_prod;
    const double* const sp_pres;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_WALL_PRODUCTION_TERMS_H
