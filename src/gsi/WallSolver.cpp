#include "WallSolver.h"

namespace Mutation{
    namespace gsi{

//-------------------------------------------------------------------------------------------

 WallSolver::WallSolver( Mutation::Thermodynamics::Thermodynamics& thermo,
                         Mutation::Transport::Transport& transport,
                         WallState* wallstate, // Pass it by reference!
                         CatalysisRateManager* catalysisrates // Pass it by reference!
                         )
                        : m_ns(thermo.nSpecies()),
                          m_transport(transport),
                          m_thermo(thermo),
                          mp_wall_state(wallstate),
                          mp_catalysis_rates(catalysisrates),
                          v_mole_fractions_wall(thermo.nSpecies(), 0.E0),
                          v_residual_function(thermo.nSpecies(), 0.E0),
                          v_unperturbed_residual_function(thermo.nSpecies(), 0.E0),
                          v_jacobian(thermo.nSpecies(), thermo.nSpecies(), 0.E0),
                          v_rhoi_wall(thermo.nSpecies(), 0.E0),
                          v_mole_fractions_edge(thermo.nSpecies(), 0.E0),
                          v_mole_fraction_gradients(thermo.nSpecies(), 0.E0),
                          v_delta_mole_fractions(thermo.nSpecies(), 0.0),
                          v_diffusion_velocities(thermo.nSpecies(), 0.E0),
                          v_wall_production_rates(thermo.nSpecies(), 0.E0),
                          m_perturbation(1.E-5),
                          zero_electric_conductivity(0.E0){

    setMaxIterations(100);
    setWriteConvergenceHistory(false);
 
}

//-------------------------------------------------------------------------------------------

void WallSolver::updateFunction( Mutation::Numerics::RealVector& lv_mole_fractions_wall ){

    for( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_mole_fractions_wall(i_ns) = lv_mole_fractions_wall(i_ns);
    }

    computePartialDensitiesfromMoleFractions( lv_mole_fractions_wall, v_rhoi_wall );
    m_thermo.setState( &v_rhoi_wall(0), &m_Twall, set_state_with_rhoi_Twall );
    mp_wall_state->setWallState( &v_rhoi_wall(0), &m_Twall );

    computeMoleFractionGradientatWall(lv_mole_fractions_wall);
    m_transport.stefanMaxwell(&v_mole_fraction_gradients(0), &v_diffusion_velocities(0), zero_electric_conductivity); 

    mp_catalysis_rates->computeRate(*mp_wall_state, &v_wall_production_rates(0));

    for(int i_ns = 0 ; i_ns < m_ns ; i_ns++){
        v_residual_function(i_ns) = v_rhoi_wall(i_ns) * v_diffusion_velocities(i_ns) - v_wall_production_rates(i_ns);
    }

}

//-------------------------------------------------------------------------------------------

void WallSolver::updateJacobian(Mutation::Numerics::RealVector& lv_mole_fractions_wall){

    double l_unperturbed_mole_fraction;
    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_unperturbed_residual_function(i_ns) = v_residual_function(i_ns);
    }

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++){

        l_unperturbed_mole_fraction = v_mole_fractions_wall(i_ns);
        lv_mole_fractions_wall(i_ns) += m_perturbation;

        updateFunction( lv_mole_fractions_wall );

        // Update Jacobian column
        for( int i_equation = 0 ; i_equation < m_ns ; i_equation++ ){
            v_jacobian(i_equation, i_ns) = (v_residual_function(i_equation) - v_unperturbed_residual_function(i_equation)) / m_perturbation ;
        }

        lv_mole_fractions_wall(i_ns) = l_unperturbed_mole_fraction;
    }

}

//-------------------------------------------------------------------------------------------

Mutation::Numerics::RealVector& WallSolver::systemSolution(){

    static Mutation::Numerics::QRP<double> m_qr(v_jacobian);
    m_qr.solve(v_delta_mole_fractions, v_unperturbed_residual_function); // Unperturbed?
    return v_delta_mole_fractions;

}

//-------------------------------------------------------------------------------------------

double WallSolver::norm(){

    return v_residual_function.normInf();

}

//-------------------------------------------------------------------------------------------

void WallSolver::setDataSolver(const double * const l_rhoi_edge, const double * const l_Twall, const double& l_dx_gradient_distance_wall_edge ){
    
    m_dx_gradient_distance_wall_edge = l_dx_gradient_distance_wall_edge;
    m_Twall = l_Twall[0];

    computeMoleFractionsfromPartialDensities(l_rhoi_edge, v_mole_fractions_edge);

}

//-------------------------------------------------------------------------------------------

void WallSolver::computeMoleFractionGradientatWall( const Mutation::Numerics::RealVector& l_mole_fractions_wall){

    for(int i_ns = 0 ; i_ns < m_ns ; i_ns++){
        v_mole_fraction_gradients(i_ns) = (l_mole_fractions_wall(i_ns) - v_mole_fractions_edge(i_ns) ) / m_dx_gradient_distance_wall_edge;
    }

}

//-------------------------------------------------------------------------------------------

void WallSolver::saveUnperturbedPressure( Mutation::Numerics::RealVector& lv_rhoi_wall ){

    m_thermo.setState( &lv_rhoi_wall(0), &m_Twall, set_state_with_rhoi_Twall);
    m_Pwall = m_thermo.P();

}

//-------------------------------------------------------------------------------------------

void WallSolver::computeMoleFractionsfromPartialDensities( Mutation::Numerics::RealVector& lv_rhoi_wall, Mutation::Numerics::RealVector& lv_mole_fractions ){

    m_thermo.setState( &lv_rhoi_wall(0), &m_Twall, set_state_with_rhoi_Twall);
    for (int i_ns = 0; i_ns < m_ns; i_ns++){
        lv_mole_fractions(i_ns) = m_thermo.X()[i_ns];
    }

}

//-------------------------------------------------------------------------------------------

void WallSolver::computeMoleFractionsfromPartialDensities(const double * const p_rhoi, Mutation::Numerics::RealVector& lv_mole_fractions ){

    m_thermo.setState( p_rhoi, &m_Twall, set_state_with_rhoi_Twall);
    for (int i_ns = 0; i_ns < m_ns; i_ns++){
        lv_mole_fractions(i_ns) = m_thermo.X()[i_ns];
    }

}

//-------------------------------------------------------------------------------------------

void WallSolver::computePartialDensitiesfromMoleFractions( Mutation::Numerics::RealVector& lv_mole_fractions, Mutation::Numerics::RealVector& lv_rhoi_wall ){

    for( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        lv_rhoi_wall(i_ns) = lv_mole_fractions(i_ns) * m_Pwall * m_thermo.speciesMw(i_ns) / ( m_Twall * Mutation::RU );
    }

}

//-------------------------------------------------------------------------------------------

void WallSolver::initializeSolutionwithPartialDensities( Mutation::Numerics::RealVector& lv_rhoi_wall ){

    saveUnperturbedPressure( lv_rhoi_wall );
    computeMoleFractionsfromPartialDensities( lv_rhoi_wall, v_mole_fractions_wall );

}

//-------------------------------------------------------------------------------------------

void WallSolver::retrieveInitialMolarFractions( Mutation::Numerics::RealVector& l_mole_fractions_wall ){

    for( int i_ns = 0 ; i_ns < m_ns ; i_ns++){
        l_mole_fractions_wall(i_ns) = v_mole_fractions_wall(i_ns);
    }

}

//-------------------------------------------------------------------------------------------

void WallSolver::returnSolutioninPartialDensities( Mutation::Numerics::RealVector& lv_solution_rhoi_wall ){

    computePartialDensitiesfromMoleFractions( v_mole_fractions_wall, lv_solution_rhoi_wall );

}

    } // namespace gsi
} // namespace Mutation

