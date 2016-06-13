#ifndef DATA_SURFACE_BALANCE_SOLVER_H
#define DATA_SURFACE_BALANCE_SOLVER_H

#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "SurfaceProperties.h"
#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

/**
 * Structure which stores the necessary inputs for the SurfaceBalanceSolver class.
 */
struct DataSurfaceBalanceSolver {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism; 
    const Mutation::Utilities::IO::XmlElement& s_node_diff_model;
    const Mutation::Utilities::IO::XmlElement& s_node_prod_terms;
    SurfaceProperties& s_surf_props;
    WallState& s_wall_state;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATA_SURFACE_BALANCE_SOLVER_H
