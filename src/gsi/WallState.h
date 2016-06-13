#ifndef WALL_STATE_H
#define WALL_STATE_H

#include "Thermodynamics.h"
#include <Eigen/Dense>

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallState{

public:
    WallState(
        const Mutation::Thermodynamics::Thermodynamics& l_thermo,
        const SurfaceProperties& l_surf_props );

    /**
     * Destructor
     */
    ~WallState();
    
    void setWallState(const double* const p_mass, const double* const p_energy, const int state_variable){
    	switch(state_variable){
    	case 0:
            setWallP(*p_mass);
            setWallT(p_energy);
            break;
    	case 1:
            setWallRhoi(p_mass);
            setWallT(p_energy);
            break;
    	default:
    		// @BD giveError();
    		exit(1);
    	}
        wallStateSet();
    }

    void getWallState(double* const p_rhoi, double* const p_rhoie, const int state_variable){
        for (int i_sp = 0; i_sp < m_ns; ++i_sp) {p_rhoi[i_sp] = getWallRhoi()(i_sp);}
        for (int i_T = 0; i_T < m_nT ; ++i_T) {p_rhoie[i_T] = getWallT()(i_T); }
    }

    void setWallRhoi(const double* const p_rhoi);
    void setWallT(const double* const p_T);
    void setWallP(const double& l_p);

    const Eigen::VectorXd& getWallRhoi() const { return v_rhoi; } // @todo const & FIX
    const Eigen::VectorXd& getWallT() const { return v_T; }
    double getWallP() const { return m_p; }

    void wallStateSet();
    bool isWallStateSet() const ;

    void getConcWallSurfState(Eigen::VectorXd& lv_wall_state) const;
//    void setSurfPropState( const Eigen::VectorXd& lv_surf_props_state ){ v_surf_props_state = lv_surf_props_state; };
//    Eigen::VectorXd getSurfPropState(){ return v_surf_props_state; }

private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceProperties& m_surf_props;

    void initializeSurfState();

    const int m_ns;
    const int m_nT;
    const int m_ns_surf;

    const int m_set_state_rhoi_T;

    Eigen::VectorXd v_rhoi;
    Eigen::VectorXd v_T;
    double m_p;

    Eigen::VectorXd v_surf_props_state;

    bool m_wall_state_set;

}; // class WallState

//======================================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // WALL_STATE_H
