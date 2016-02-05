#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

WallState::WallState( const Mutation::Thermodynamics::Thermodynamics& l_thermo,
                      const SurfaceProperties& l_surf_props )
                     : m_thermo ( l_thermo ),
                       m_surf_props( l_surf_props ),
                       m_ns( l_thermo.nSpecies() ),
                       m_nT( l_thermo.nEnergyEqns() ),
                       m_ns_surf( l_surf_props.nSpeciesWall() ),
                       m_set_state_rhoi_T( 1 ),
                       v_rhoi( m_ns ),
                       v_T( m_nT ),
                       m_wall_state_set( 0 ),
                       v_surf_props_state( m_ns_surf ) {

	initializeSurfState();

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

void WallState::getConcWallSurfState( Eigen::VectorXd& lv_wall_state ) const {
	assert(lv_wall_state.size() == m_ns+ m_ns_surf );

	m_thermo.convert<Mutation::Thermodynamics::RHO_TO_CONC>( &v_rhoi[0], &lv_wall_state[0] );
	lv_wall_state.head( m_ns ) *= Mutation::NA ;
	lv_wall_state.tail( m_ns_surf ) = v_surf_props_state.head( m_ns_surf );
}

//======================================================================================

void WallState::initializeSurfState(){
	int n_sites = m_surf_props.nSites();
	double ln_total_sites = m_surf_props.nTotalSites();

	int n_sp_in_site = 0;
	double n_frac_site = 0.0;
    double n_sites_dens = 0.0;
    int pos_in_surf_props = 0;

    for ( int i_sites = 0; i_sites < n_sites; ++i_sites ){
    	n_frac_site = m_surf_props.fracSite( i_sites );
    	n_sp_in_site = m_surf_props.nSpeciesSite( i_sites );
    	n_sites_dens = ln_total_sites * n_frac_site / n_sp_in_site;
        for ( int i_sp_sites = 0; i_sp_sites < n_sp_in_site; i_sp_sites++ ){
        	v_surf_props_state[ pos_in_surf_props ] = n_sites_dens;
        	pos_in_surf_props++;
        }
	}

}

//======================================================================================

    } // GasSurfaceInteraction
} // Mutation
