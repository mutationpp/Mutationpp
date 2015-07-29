#ifndef SURFACEDESCRIPTION_H
#define SURFACEDESCRIPTION_H

#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"
#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceDescription{

public:
    SurfaceDescription( const Mutation::Thermodynamics::Thermodynamics& l_thermo, const std::string& l_gsi_mechanism, const Mutation::Utilities::IO::XmlElement& l_node_surf_props );
    ~SurfaceDescription();

    void setWallState( const double* const p_rhoi, const double* const p_rhoie, const int state_variable );
    bool isWallStateSet() const;

    const Mutation::Numerics::RealVector& getWallRhoi() const { return mp_wall_state->getWallRhoi(); }
    const Mutation::Numerics::RealVector& getWallT() const { return mp_wall_state->getWallT(); }

private:
    WallState* mp_wall_state;
    SurfaceProperties* mp_surf_props;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACEDESCRIPTION_H
