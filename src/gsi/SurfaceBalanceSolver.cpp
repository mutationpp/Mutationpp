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
                                            v_mole_frac_wall( m_thermo.nSpecies(), 0.E0 ),
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

//    std::fill( &v_sep_mass_prod_rate.begin(), &v_sep_mass_prod_rate.end(), 0.0 ); /
    for ( int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){
        v_sep_mass_prod_rate( i_ns ) = 0.E0;
    }

    for ( int i_prod_terms = 0 ; i_prod_terms < v_surf_prod.size() ; ++i_prod_terms ){
        v_surf_prod[i_prod_terms]->productionRate( v_sep_mass_prod_rate );
    }

    for ( int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){
        lv_mass_prod_rate( i_ns ) += v_sep_mass_prod_rate( i_ns );
    }

}

//======================================================================================

void SurfaceBalanceSolver::setDiffusionModel( const Mutation::Numerics::RealVector& lv_mole_frac_edge, const double& l_dx ){

    mp_diff_vel_calc->setDiffusionModel( lv_mole_frac_edge, l_dx );

}

//======================================================================================

void SurfaceBalanceSolver::solveSurfaceBalance(){

    // errorUninitializedDiffusionModel

    // errorWallStateNotSetYet
    
    // Get Wall State

    //    solve();

}

//======================================================================================

void SurfaceBalanceSolver::updateFunction( Mutation::Numerics::RealVector& lv_mole_frac_wall ) {

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ) {
        v_mole_frac_wall( i_ns ) = lv_mole_frac_wall( i_ns );
    }

    // computePartialDensitiesfromMoleFractions
    // m_thermo.setState
    // set_wall_state
    //
    // gradient
    // stefanMaxwell
    //
    // loop over v_surf_prod->computeRate(v_wall_prod_rates)
    //
    // v_residual_function


}

//======================================================================================

void SurfaceBalanceSolver::updateJacobian( Mutation::Numerics::RealVector& lv_mole_frac_wall ) {

    // Backup unperturbed mole fraction
    // Backup Unperturbed residual function
    
    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++){

        // Perturb the mole fractions
        // l_unperturbed_mole_fraction = v_mole_fractions_wall(i_ns);
        // lv_mole_fractions_wall(i_ns) += m_perturbation;

        // updateFunction( lv_mole_fractions_wall );

        // Update Jacobian column
//        for( int i_equation = 0 ; i_equation < m_ns ; i_equation++ ){
//            v_jacobian(i_equation, i_ns) = ( v_residual_function( i_equation ) - v_unperturbed_residual_function( i_equation ) ) / m_perturbation ;
//        }

        // Unperturbe mole fractions
        // lv_mole_fractions_wall(i_ns) = l_unperturbed_mole_fraction;
        // Or equivalently lv_mole_fractions_wall(i_ns) = v_mole_fractions_wall(i_ns);

    }


}

//======================================================================================

Mutation::Numerics::RealVector& SurfaceBalanceSolver::systemSolution(){

//    static Mutation::Numerics::QRP<double> m_qr(v_jacobian);
//    m_qr.solve(v_d_mole_frac, v_unperturbed_residual_function); // Unperturbed?
    return v_d_mole_frac;

}

//======================================================================================

double SurfaceBalanceSolver::norm(){

    return 0.0;//v_residual_function.normInf();

}

//======================================================================================

void SurfaceBalanceSolver::addSurfaceProductionTerm( WallProductionTerms* lp_wall_prod_term ){

    v_surf_prod.push_back( lp_wall_prod_term );

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
