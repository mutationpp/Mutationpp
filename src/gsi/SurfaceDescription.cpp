#include <stddef.h>

#include "Utilities.h"

#include "SurfaceDescription.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

SurfaceDescription::SurfaceDescription( const Mutation::Thermodynamics::Thermodynamics& l_thermo, const std::string& l_gsi_mechanism, const Mutation::Utilities::IO::XmlElement& l_node_surf_props ) 
                                      : mp_surf_props( NULL ),
                                        mp_wall_state( NULL ) {
 
    mp_surf_props = Mutation::Utilities::Config::Factory<SurfaceProperties>::create( l_gsi_mechanism, l_node_surf_props );
    mp_wall_state = new WallState( l_thermo /*, mp_surf_props*/ );

}

//======================================================================================

SurfaceDescription::~SurfaceDescription() {

    if ( mp_wall_state != NULL ) { delete mp_wall_state; }
    if ( mp_surf_props != NULL ) { delete mp_surf_props; }

}

//======================================================================================

void SurfaceDescription::setWallState( const double* const p_rhoi, const double* const p_rhoie, const int state_variable ){
    
    mp_wall_state->setWallRhoi( p_rhoi );
    mp_wall_state->setWallT( p_rhoie );

    mp_wall_state->wallStateSet();

}

//======================================================================================

void SurfaceDescription::getWallState( double* const p_rhoi, double* const p_rhoie, const int state_variable ){

    for ( int i_sp = 0; i_sp < mp_wall_state->getWallRhoi().size() ; ++i_sp ) { /** @todo This might be really really slow... */
        p_rhoi[ i_sp ] = mp_wall_state->getWallRhoi()( i_sp );
    }

    for ( int i_T = 0; i_T < mp_wall_state->getWallRhoi().size() ; ++i_T ) {
        p_rhoie[ i_T ] = mp_wall_state->getWallT()( i_T );
    }

}

//======================================================================================

bool SurfaceDescription::isWallStateSet() const {

    return mp_wall_state->isWallStateSet();

}

//======================================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
