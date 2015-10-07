#ifndef WALLSTATE_H
#define WALLSTATE_H

#include "Thermodynamics.h"
#include <Eigen/Dense>

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallState{

public:
    WallState( const Mutation::Thermodynamics::Thermodynamics& thermo );
    ~WallState();
    
    void setWallRhoi( const double* const p_rhoi );
    void setWallT( const double* const p_T );

    const Eigen::VectorXd& getWallRhoi() const { return v_rhoi; }
    const Eigen::VectorXd& getWallT() const { return v_T; }

    void wallStateSet();
    bool isWallStateSet() const ;

private:
    const int m_ns;
    const int m_nT;

    Eigen::VectorXd v_rhoi;
    Eigen::VectorXd v_T;

    bool m_wall_state_set;

};

//======================================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // WALLSTATE_H
