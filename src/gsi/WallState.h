#ifndef WALLSTATE_H
#define WALLSTATE_H

#include "Thermodynamics.h"
#include <Eigen/Dense>

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallState{

public:
    WallState( const Mutation::Thermodynamics::Thermodynamics& l_thermo,
               const SurfaceProperties& l_surf_props );
    ~WallState();
    
    void setWallState( const double* const p_rhoi, const double* const p_rhoie, const int state_variable ){
        setWallRhoi( p_rhoi );
        setWallT( p_rhoie );
        wallStateSet();
    }

    void getWallState( double* const p_rhoi, double* const p_rhoie, const int state_variable ){
        for ( int i_sp = 0; i_sp < m_ns ; ++i_sp ) { p_rhoi[ i_sp ] = getWallRhoi()( i_sp ); }
        for ( int i_T = 0; i_T < m_nT ; ++i_T ) { p_rhoie[ i_T ] = getWallT()( i_T ); }
    }

    void setWallRhoi( const double* const p_rhoi );
    void setWallT( const double* const p_T );

    const Eigen::VectorXd& getWallRhoi() const { return v_rhoi; }
    const Eigen::VectorXd& getWallT() const { return v_T; }

    void wallStateSet();
    bool isWallStateSet() const ;

    void setSurfPropState( const Eigen::VectorXd& lv_surf_props_state ){ v_surf_props_state = lv_surf_props_state; };
    Eigen::VectorXd getSurfPropState(){ return v_surf_props_state; }

private:
    const int m_ns;
    const int m_nT;

    Eigen::VectorXd v_rhoi;
    Eigen::VectorXd v_T;

    Eigen::VectorXd v_surf_props_state;

    bool m_wall_state_set;

};

//======================================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // WALLSTATE_H
