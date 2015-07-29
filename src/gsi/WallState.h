#ifndef WALLSTATE_H
#define WALLSTATE_H

#include "Numerics.h"
#include "Thermodynamics.h"
#include "Vector.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallState{

public:
    WallState( const Mutation::Thermodynamics::Thermodynamics& thermo );
    ~WallState();
    
    void setWallRhoi( const double* const p_rhoi );
    void setWallT( const double* const p_T );

    const Mutation::Numerics::RealVector& getWallRhoi() const { return v_rhoi; }
    const Mutation::Numerics::RealVector& getWallT() const { return v_T; }

    void wallStateSet();
    bool isWallStateSet() const ;

private:
    const int m_ns;
    const int m_nT;

    Mutation::Numerics::RealVector v_rhoi;
    Mutation::Numerics::RealVector v_T;

    bool m_wall_state_set;

};

//======================================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // WALLSTATE_H
