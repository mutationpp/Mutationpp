#include "Errors.h"
#include "Thermodynamics.h"

#include "SurfaceProperties.h"
#include "WallState.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//==============================================================================

WallState::WallState(
    const Mutation::Thermodynamics::Thermodynamics& l_thermo,
    const SurfaceProperties& l_surf_props)
    : m_thermo(l_thermo),
      m_surf_props(l_surf_props),
      m_ns(l_thermo.nSpecies()),
      m_nT(l_thermo.nEnergyEqns()),
      m_ns_surf(l_surf_props.nWallSpecies()),
      m_set_state_rhoi_T(1),
      mv_rhoi(m_ns),
      mv_T(m_nT),
      m_is_wall_state_set(false),
      mv_surf_props_state(m_ns_surf)
{
	initializeSurfState();
}

//==============================================================================

WallState::~WallState(){ }

//==============================================================================

void WallState::setWallState(
    const double* const p_mass,
    const double* const p_energy,
    const int state_var)
{
  	switch(state_var){
    case 0:
        setWallP(*p_mass);
        setWallT(p_energy);
        break;
    case 1:
        setWallRhoi(p_mass);
        setWallT(p_energy);
        break;
    default:
        throw InvalidInputError("variable set", state_var)
        << "This variable-set is not implemented in setWallState"
        << ". Possible variable-sets are:\n"
        << "  0: (pressure, temperature)\n"
        << "  1: (species densities, temperature)\n";
    }
    m_is_wall_state_set = true;
}

//==============================================================================

void WallState::getWallState(
    double* const p_rhoi,
    double* const p_rhoie,
    const int state_var)
{
  	switch(state_var){
    case 1:
        for (int i_sp = 0; i_sp < m_ns; ++i_sp){
            p_rhoi[i_sp] = getWallRhoi()(i_sp);
        }
        for (int i_T = 0; i_T < m_nT ; ++i_T) {
            p_rhoie[i_T] = getWallT()(i_T);
        }
        break;
    default:
        throw InvalidInputError("variable get", state_var)
        << "This variable-get is not implemented in getWallState"
        << ". Possible variable-sets are:\n"
        << "  1: (species densities, temperature)\n";
  	}
}

//==============================================================================

void WallState::setWallRhoi(const double* const p_rhoi){
	mv_rhoi = Eigen::Map<const Eigen::VectorXd>(p_rhoi, m_ns);
}

//==============================================================================

void WallState::setWallT(const double* const p_T){
    mv_T = Eigen::Map<const Eigen::VectorXd>(p_T, m_nT);
}

//==============================================================================

void WallState::setWallP(const double& p){
    m_p = p;
}

//==============================================================================

void WallState::getNdStateGasSurf(Eigen::VectorXd& v_wall_state) const
{
	assert(v_wall_state.size() == m_ns+ m_ns_surf);

	m_thermo.convert<Mutation::Thermodynamics::RHO_TO_CONC>(
        mv_rhoi.data(),
        v_wall_state.data());

	v_wall_state.head(m_ns) *= Mutation::NA ;
	v_wall_state.tail(m_ns_surf) = mv_surf_props_state.head(m_ns_surf);

}

//==============================================================================

void WallState::initializeSurfState()
{
	size_t n_sites = m_surf_props.nSiteCategories();
	double n_total_sites = m_surf_props.nTotalSites();

	int n_sp_in_site = 0;
	double n_frac_site = 0.0;
    double n_sites_dens = 0.0;
    int pos_in_surf_props = 0;

    for (int i_sites = 0; i_sites < n_sites; ++i_sites){
    	n_frac_site = m_surf_props.fracSite(i_sites);
    	n_sp_in_site = m_surf_props.nSpeciesinSite(i_sites);
    	n_sites_dens = n_total_sites * n_frac_site / n_sp_in_site;
        for (int i_sp_sites = 0; i_sp_sites < n_sp_in_site; i_sp_sites++){
        	mv_surf_props_state[pos_in_surf_props] = n_sites_dens;
        	pos_in_surf_props++;
        }
	}
}

//==============================================================================

    } // GasSurfaceInteraction
} // Mutation
