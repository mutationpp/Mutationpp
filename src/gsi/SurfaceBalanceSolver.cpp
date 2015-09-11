#include "DataWallProductionTerms.h"
#include "SurfaceBalanceSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

SurfaceBalanceSolver::SurfaceBalanceSolver( Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport, const std::string& l_gsi_mechanism, const Mutation::Utilities::IO::XmlElement& l_node_diff_model, const Mutation::Utilities::IO::XmlElement& l_node_prod_terms, SurfaceDescription& l_surf_descr )
                                          : m_thermo( l_thermo ),
                                            m_surf_descr( l_surf_descr ),
                                            mp_diff_vel_calc( NULL ),
                                            mp_mass_blowing_rate( NULL ),
                                            m_ns( m_thermo.nSpecies() ),
                                            v_rhoi( m_thermo.nSpecies(), 0.E0 ),
                                            v_work( m_thermo.nSpecies(), 0.E0 ),
                                            v_X( m_thermo.nSpecies(), 0.E0 ),
                                            v_dX( m_thermo.nSpecies(), 0.E0 ),
                                            v_f( m_thermo.nSpecies(), 0.E0 ),
                                            v_f_unpert( m_thermo.nSpecies(), 0.E0 ),
                                            v_jac(m_thermo.nSpecies(), m_thermo.nSpecies(), 0.E0),
                                            m_pert( 1.E-5 ),
                                            set_state_with_rhoi_T( 1 ),
                                            v_sep_mass_prod_rate( m_thermo.nSpecies(), 0.E0 ) {
                                            

    Mutation::Utilities::IO::XmlElement::const_iterator iter_prod_terms = l_node_prod_terms.begin();

//    { /** @todo Make me a function */
    DataWallProductionTerms l_data_wall_production_terms = { m_thermo, l_gsi_mechanism, *iter_prod_terms, m_surf_descr }; 
    for( ; iter_prod_terms != l_node_prod_terms.end() ; ++iter_prod_terms ){
        addSurfaceProductionTerm( Mutation::Utilities::Config::Factory<WallProductionTerms>::create( iter_prod_terms->tag(), l_data_wall_production_terms ) );
    }
    errorEmptyWallProductionTerms();
//    } // Make me a function
 
    // DiffusionVelocityCalculator
    mp_diff_vel_calc = new DiffusionVelocityCalculator( m_thermo, l_transport );

    // MassBlowingRate
    mp_mass_blowing_rate = new MassBlowingRate();

    // Setup NewtonSolver
    setMaxIterations(3);
    setWriteConvergenceHistory(false);

}

//======================================================================================

SurfaceBalanceSolver::~SurfaceBalanceSolver(){

    for ( std::vector<WallProductionTerms*>::iterator iter = v_surf_prod.begin() ; iter != v_surf_prod.end() ; ++iter ){
        delete (*iter);
    }
    v_surf_prod.clear();
    if ( mp_diff_vel_calc != NULL ){ delete mp_diff_vel_calc; }
    if ( mp_mass_blowing_rate != NULL ){ delete mp_mass_blowing_rate; }

}

//======================================================================================

void SurfaceBalanceSolver::setWallState( const double* const l_mass, const double* const l_energy, const int state_variable ){
    m_surf_descr.setWallState( l_mass, l_energy, state_variable );
}

//======================================================================================

void SurfaceBalanceSolver::computeGSIProductionRate( Mutation::Numerics::RealVector& lv_mass_prod_rate ){

    errorWallStateNotSet();
    for ( int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){ v_sep_mass_prod_rate( i_ns ) = 0.E0; }

    for ( int i_prod_terms = 0 ; i_prod_terms < v_surf_prod.size() ; ++i_prod_terms ){
        v_surf_prod[i_prod_terms]->productionRate( v_sep_mass_prod_rate );
    }

    for ( int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){ /////// WTF????????????
        lv_mass_prod_rate( i_ns ) = v_sep_mass_prod_rate( i_ns );
    }

}

//======================================================================================

void SurfaceBalanceSolver::setDiffusionModel( const Mutation::Numerics::RealVector& lv_mole_frac_edge, const double& l_dx ){

    mp_diff_vel_calc->setDiffusionModel( lv_mole_frac_edge, l_dx );

}

//======================================================================================

void SurfaceBalanceSolver::solveSurfaceBalance( const Mutation::Numerics::RealVector& lv_rhoi, const Mutation::Numerics::RealVector& lv_T ){

    // errorUninitializedDiffusionModel
    // errorWallStateNotSetYet

    v_rhoi = lv_rhoi;
    m_Twall = lv_T(0);
    saveUnperturbedPressure( v_rhoi );

    computeMoleFracfromPartialDens( v_rhoi, v_X );
    v_X = solve( v_X );
    computePartialDensfromMoleFrac( v_X, v_rhoi );

    setWallState( &v_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    
}

//======================================================================================

void SurfaceBalanceSolver::updateFunction( Mutation::Numerics::RealVector& lv_mole_frac ) {

    computePartialDensfromMoleFrac( lv_mole_frac, v_rhoi );
    m_thermo.setState( &v_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    setWallState( &v_rhoi(0), &m_Twall, set_state_with_rhoi_T );

    mp_diff_vel_calc->computeDiffusionVelocities( lv_mole_frac, v_work );

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_f(i_ns) = v_rhoi(i_ns) * v_work(i_ns); 
    }

    computeGSIProductionRate( v_work );

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_f(i_ns) -= v_work(i_ns); 
    }

}

//======================================================================================

void SurfaceBalanceSolver::updateJacobian( Mutation::Numerics::RealVector& lv_mole_frac ) {

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_f_unpert(i_ns) = v_f(i_ns);
    }
    
    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++){

        m_X_unpert = lv_mole_frac(i_ns);
        lv_mole_frac(i_ns) += m_pert;

        updateFunction( lv_mole_frac );

        // Update Jacobian column
        for( int i_eq = 0 ; i_eq < m_ns ; i_eq++ ){
            v_jac(i_eq, i_ns) = ( v_f( i_eq ) - v_f_unpert( i_eq ) ) / m_pert ;
        }

        // Unperturbe mole fractions
        lv_mole_frac(i_ns) = m_X_unpert;

    }

}

//======================================================================================

Mutation::Numerics::RealVector& SurfaceBalanceSolver::systemSolution(){

    static Mutation::Numerics::QRP<double> m_qr(v_jac);
    m_qr.solve(v_dX, v_f_unpert); // Unperturbed?
    return v_dX;

}

//======================================================================================

double SurfaceBalanceSolver::norm(){

    return v_f.normInf();

}

//======================================================================================

void SurfaceBalanceSolver::addSurfaceProductionTerm( WallProductionTerms* lp_wall_prod_term ){

    v_surf_prod.push_back( lp_wall_prod_term );

}

//======================================================================================

void SurfaceBalanceSolver::saveUnperturbedPressure( Mutation::Numerics::RealVector& lv_rhoi ){

    m_thermo.setState( &lv_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    m_Pwall = m_thermo.P();

}

//======================================================================================

void SurfaceBalanceSolver::computeMoleFracfromPartialDens( Mutation::Numerics::RealVector& lv_rhoi, Mutation::Numerics::RealVector& lv_xi ){

    m_thermo.setState( &lv_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        lv_xi(i_ns) = m_thermo.X()[i_ns];
    }

}

//======================================================================================

void SurfaceBalanceSolver::computePartialDensfromMoleFrac( Mutation::Numerics::RealVector& lv_xi, Mutation::Numerics::RealVector& lv_rhoi ){

    for( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        lv_rhoi(i_ns) = lv_xi(i_ns) * m_Pwall * m_thermo.speciesMw(i_ns) / ( m_Twall * Mutation::RU );
    }

}
//======================================================================================

void SurfaceBalanceSolver::errorEmptyWallProductionTerms() const {

    if ( v_surf_prod.size() == 0 ){

        std::cerr << "At least one wall production term should be provided in the input file!" << std::endl;
        exit(1);

    }

}

//======================================================================================

void SurfaceBalanceSolver::errorWallStateNotSet() const {
    
    if ( m_surf_descr.isWallStateSet() == 0 ){ 
        std::cerr << "The wall state must have been set!" << std::endl; /** @todo better error */
        exit(1);
    }

}

//======================================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
