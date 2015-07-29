#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

WallState::WallState( const Mutation::Thermodynamics::Thermodynamics& thermo )
                     : m_ns( thermo.nSpecies() ),
                       m_nT( thermo.nEnergyEqns() ),
                       v_rhoi( m_ns),
                       v_T( m_nT ),
                       m_wall_state_set( 0 ) {

}

//======================================================================================

WallState::~WallState(){ }

//======================================================================================

void WallState::setWallRhoi( const double* const p_rhoi ){

    for (int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_rhoi( i_ns ) = p_rhoi[ i_ns ];
    }

}

//======================================================================================

void WallState::setWallT( const double* const p_T ){

    for (int i_nT = 0 ; i_nT < m_nT ; i_nT++ ){
        v_T( i_nT ) = p_T[ i_nT ];
    }

}

//======================================================================================

void WallState::wallStateSet(){
    m_wall_state_set = 1;
}

//======================================================================================

bool WallState::isWallStateSet() const {
    return m_wall_state_set;
}

//======================================================================================

    } // GasSurfaceInteraction
} // Mutation
