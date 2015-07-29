#include <stddef.h>

#include "DiffusionVelocityCalculator.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

DiffusionVelocityCalculator::DiffusionVelocityCalculator( const Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport )
                                                        : v_rhs( l_thermo.nSpecies(), 0.0),
                                                          m_transport( l_transport ) {

    v_driving_forces.push_back( new DiffusionDrivingForces( l_thermo ) );

}

//======================================================================================

DiffusionVelocityCalculator::~DiffusionVelocityCalculator(){

    delete v_driving_forces[0];

}

//======================================================================================

void DiffusionVelocityCalculator::setDiffusionModel( const Mutation::Numerics::RealVector& lv_mole_frac_edge, const double& l_dx ){
    v_driving_forces[0]->setDiffusionCalculator( lv_mole_frac_edge, l_dx);
}

//======================================================================================

void DiffusionVelocityCalculator::computeDiffusionVelocities( const Mutation::Numerics::RealVector& lv_mole_frac, Mutation::Numerics::RealVector& lv_diff_velocities ){

    v_rhs = 0.E0;

    for ( int i_driv = 0 ; i_driv < v_driving_forces.size() ; ++i_driv ){
        v_driving_forces[i_driv]->computeDrivingForces( lv_mole_frac, v_rhs );
    }

    double electric_field = 0.E0;

    m_transport.stefanMaxwell( &v_rhs(0), &lv_diff_velocities(0), electric_field ); 

}

//======================================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
