#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

WallState::WallState(
    const Mutation::Thermodynamics::Thermodynamics& l_thermo,
    const SurfaceProperties& l_surf_props)
    : m_thermo(l_thermo),
      m_surf_props(l_surf_props),
      m_ns(l_thermo.nSpecies()),
      m_nT(l_thermo.nEnergyEqns()),
      m_ns_surf(l_surf_props.nWallSpecies()),
      m_set_state_rhoi_T(1),
      v_rhoi(m_ns),
      v_T(m_nT),
      m_wall_state_set(false),
      v_surf_props_state(m_ns_surf)
{
	initializeSurfState();
}

//======================================================================================

WallState::~WallState(){ }

//======================================================================================

void WallState::setWallRhoi(const double* const p_rhoi){
	v_rhoi = Eigen::Map<const Eigen::VectorXd>(p_rhoi, m_ns);
}

//======================================================================================

void WallState::setWallT(const double* const p_T){
    v_T = Eigen::Map<const Eigen::VectorXd>(p_T, m_nT);
}

//======================================================================================

void WallState::setWallP(const double& l_p){ //@todo Pass it by pointer
    m_p = l_p;
}

//======================================================================================

void WallState::wallStateSet(){
    m_wall_state_set = true;
}

//======================================================================================

bool WallState::isWallStateSet() const {
    return m_wall_state_set;
}

//======================================================================================

void WallState::getNdStateGasSurf(Eigen::VectorXd& lv_wall_state) const
{
	assert(lv_wall_state.size() == m_ns+ m_ns_surf);

	m_thermo.convert<Mutation::Thermodynamics::RHO_TO_CONC>(v_rhoi.data(), lv_wall_state.data());
	lv_wall_state.head(m_ns) *= Mutation::NA ;
	lv_wall_state.tail(m_ns_surf) = v_surf_props_state.head(m_ns_surf);

}

//======================================================================================

void WallState::initializeSurfState()
{
	size_t n_sites = m_surf_props.nSiteCategories();
	double ln_total_sites = m_surf_props.nTotalSites();

	int n_sp_in_site = 0;
	double n_frac_site = 0.0;
    double n_sites_dens = 0.0;
    int pos_in_surf_props = 0;

    for (int i_sites = 0; i_sites < n_sites; ++i_sites){
    	n_frac_site = m_surf_props.fracSite(i_sites);
    	n_sp_in_site = m_surf_props.nSpeciesinSite(i_sites);
    	n_sites_dens = ln_total_sites * n_frac_site / n_sp_in_site;
        for (int i_sp_sites = 0; i_sp_sites < n_sp_in_site; i_sp_sites++){
        	v_surf_props_state[pos_in_surf_props] = n_sites_dens;
        	pos_in_surf_props++;
        }
	}
}

//======================================================================================

    } // GasSurfaceInteraction
} // Mutation
