#include "WallSteadyStateCoverageSolver.h"

namespace Mutation{
    namespace gsi{

//=============================================================================
//=============================================================================
//=============================================================================

WallSteadyStateCoverageSolver::WallSteadyStateCoverageSolver( const Mutation::Thermodynamics::Thermodynamics& l_thermo, WallState& l_wall_state )
                                                            : m_wall_state(l_wall_state),
                                                              v_residual_function( l_thermo.nSpecies(), 0.E0),
                                                              v_jacobian( l_thermo.nSpecies(), l_thermo.nSpecies(), 0.E0 )
{
}
//=============================================================================
WallSteadyStateCoverageSolver::~WallSteadyStateCoverageSolver(){ }
//=============================================================================
void WallSteadyStateCoverageSolver::updateFunction(){



// computeProductionRates();

// update v_residual_function
//    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
//        v_residual_function(i_ns) = 0.0;
//    }

}
//=============================================================================
void WallSteadyStateCoverageSolver::updateJacobian( Mutation::Numerics::RealVector& ){

// Do perturbation

// Recompute everything

// Update Jacobian

// Undo perturbation

// Exit

 }
//=============================================================================
Mutation::Numerics::RealVector& WallSteadyStateCoverageSolver::WallSteadyStateCoverageSolver::systemSolution(){
//    m_qr.solve(v_delta_mole_fractions, v_unperturbed_residual_function); // Unperturbed?
//    return v_delta_mole_fractions;
}
//=============================================================================
double WallSteadyStateCoverageSolver::WallSteadyStateCoverageSolver::norm(){ 
    return 0.0;//v_residual_function.normInf();
}

//=============================================================================
//=============================================================================
//=============================================================================

    } // namespace gsi
} // namespace Mutation
