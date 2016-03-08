#include "DataWallProductionTerms.h"
#include "SurfaceBalanceSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

class SurfaceBalanceSolverGamma : public SurfaceBalanceSolver, public Mutation::Numerics::NewtonSolver<Eigen::VectorXd, SurfaceBalanceSolverGamma> {

//======================================================================================

public:
SurfaceBalanceSolverGamma( ARGS l_data_surface_balance_solver )
                                          : m_thermo( l_data_surface_balance_solver.s_thermo ),
                                            m_surf_props( l_data_surface_balance_solver.s_surf_props ),
                                            m_wall_state( l_data_surface_balance_solver.s_wall_state ),
                                            mp_diff_vel_calc( NULL ),
                                            mp_mass_blowing_rate( NULL ),
                                            m_ns( m_thermo.nSpecies() ),
                                            v_rhoi( m_thermo.nSpecies() ), // Removed Initialization to zero
                                            v_work( m_thermo.nSpecies() ), // Removed Initialization to zero
                                            v_X( m_thermo.nSpecies() ), // Removed Initialization to zero
                                            v_dX( m_thermo.nSpecies() ), // Removed Initialization to zero
                                            v_f( m_thermo.nSpecies() ), // Removed Initialization to zero
                                            v_f_unpert( m_thermo.nSpecies() ), // Removed Initialization to zero
                                            v_jac(m_thermo.nSpecies(), m_thermo.nSpecies() ), // Removed Initialization to zero
                                            m_pert( 1.E-5 ),
                                            set_state_with_rhoi_T( 1 ),
                                            v_sep_mass_prod_rate( m_thermo.nSpecies() ) { // Removed Initialization to zero
                                            

    Mutation::Utilities::IO::XmlElement::const_iterator iter_prod_terms = l_data_surface_balance_solver.s_node_prod_terms.begin();

    DataWallProductionTerms l_data_wall_production_terms = { m_thermo, l_data_surface_balance_solver.s_gsi_mechanism, *iter_prod_terms, m_surf_props, m_wall_state }; 
    for( ; iter_prod_terms != l_data_surface_balance_solver.s_node_prod_terms.end() ; ++iter_prod_terms ){
        addSurfaceProductionTerm( Mutation::Utilities::Config::Factory<WallProductionTerms>::create( iter_prod_terms->tag(), l_data_wall_production_terms ) );
    }
    errorEmptyWallProductionTerms();
 
    // DiffusionVelocityCalculator
    mp_diff_vel_calc = new DiffusionVelocityCalculator( m_thermo, l_data_surface_balance_solver.s_transport );

    // MassBlowingRate
    DataMassBlowingRate l_data_mass_blowing_rate = { *v_surf_prod[0] };
  //  mp_mass_blowing_rate = Mutation::Utilities::Config::Factory<MassBlowingRate>::create("zero", l_data_mass_blowing_rate);
    mp_mass_blowing_rate = Mutation::Utilities::Config::Factory<MassBlowingRate>::create("ablation", l_data_mass_blowing_rate);

    // Setup NewtonSolver
    setMaxIterations(3); // Write a wrapper around the maximum number of iterations.
    setWriteConvergenceHistory(false);

}

//======================================================================================

~SurfaceBalanceSolverGamma(){

    for ( std::vector<WallProductionTerms*>::iterator iter = v_surf_prod.begin() ; iter != v_surf_prod.end() ; ++iter ){
        delete (*iter);
    }
    v_surf_prod.clear();
    if ( mp_diff_vel_calc != NULL ){ delete mp_diff_vel_calc; }
    if ( mp_mass_blowing_rate != NULL ){ delete mp_mass_blowing_rate; }

}

//======================================================================================

void computeGSIProductionRate( Eigen::VectorXd& lv_mass_prod_rate ){

    errorWallStateNotSet();
    for ( int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){ v_sep_mass_prod_rate( i_ns ) = 0.E0; } // v_sep_mass_prod_rate.setZero();

    for ( int i_prod_terms = 0 ; i_prod_terms < v_surf_prod.size() ; ++i_prod_terms ){
        v_surf_prod[i_prod_terms]->productionRate( v_sep_mass_prod_rate );
    }

    for ( int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){ // lv_mass_prod_rate = v_sep_mass_prod_rate;
        lv_mass_prod_rate( i_ns ) = v_sep_mass_prod_rate( i_ns );
    }

}

//======================================================================================

void setDiffusionModel( const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx ){
    mp_diff_vel_calc->setDiffusionModel( lv_mole_frac_edge, l_dx );
}

//======================================================================================

void solveSurfaceBalance( const Eigen::VectorXd& lv_rhoi, const Eigen::VectorXd& lv_T ){

    // errorUninitializedDiffusionModel
    // errorWallStateNotSetYet

    v_rhoi = lv_rhoi;
    m_Twall = lv_T(0);
    saveUnperturbedPressure( v_rhoi );

    computeMoleFracfromPartialDens( v_rhoi, v_X );
    v_X = solve( v_X );
    computePartialDensfromMoleFrac( v_X, v_rhoi );

    m_wall_state.setWallState( &v_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    
}

//======================================================================================

void updateFunction( Eigen::VectorXd& lv_mole_frac ) {

    computePartialDensfromMoleFrac( lv_mole_frac, v_rhoi );
    m_thermo.setState( &v_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    m_wall_state.setWallState( &v_rhoi(0), &m_Twall, set_state_with_rhoi_T );

    mp_diff_vel_calc->computeDiffusionVelocities( lv_mole_frac, v_work );

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_f(i_ns) = v_rhoi(i_ns) * v_work(i_ns); 
    }

    computeGSIProductionRate( v_work );

    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        v_f(i_ns) -= v_work(i_ns); 
    }

    //@totuesday add massblowing rate
     double l_blowing = mp_mass_blowing_rate->computeBlowingFlux() / v_rhoi.sum();

     for (int i_ns = 0 ; i_ns < m_ns ; ++i_ns ){
         v_f(i_ns) -= l_blowing * v_rhoi(i_ns);
     }

}

//======================================================================================

void updateJacobian( Eigen::VectorXd& lv_mole_frac ) {

        v_f_unpert = v_f;
    
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

Eigen::VectorXd& systemSolution(){
    v_dX = v_jac.partialPivLu().solve(v_f_unpert);
    return v_dX;
}

//======================================================================================

double norm(){
    return v_f.lpNorm<Eigen::Infinity>();
}

//======================================================================================
private:
void addSurfaceProductionTerm( WallProductionTerms* lp_wall_prod_term ){
    v_surf_prod.push_back( lp_wall_prod_term );
}

//======================================================================================

void saveUnperturbedPressure( Eigen::VectorXd& lv_rhoi ){
    m_thermo.setState( &lv_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    m_Pwall = m_thermo.P();
}

//======================================================================================

void computeMoleFracfromPartialDens( Eigen::VectorXd& lv_rhoi, Eigen::VectorXd& lv_xi ){
    m_thermo.setState( &lv_rhoi(0), &m_Twall, set_state_with_rhoi_T );
    for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        lv_xi(i_ns) = m_thermo.X()[i_ns];
    }
}

//======================================================================================

void computePartialDensfromMoleFrac( Eigen::VectorXd& lv_xi, Eigen::VectorXd& lv_rhoi ){
    for( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
        lv_rhoi(i_ns) = lv_xi(i_ns) * m_Pwall * m_thermo.speciesMw(i_ns) / ( m_Twall * Mutation::RU );
    }
}
//======================================================================================

void errorEmptyWallProductionTerms() const {
    if ( v_surf_prod.size() == 0 ){
        std::cerr << "At least one wall production term should be provided in the input file!" << std::endl;
        exit(1);
    }
}

//======================================================================================

void errorWallStateNotSet() const {
    if ( m_wall_state.isWallStateSet() == 0 ){ 
        std::cerr << "The wall state must have been set!" << std::endl; /** @todo better error */
        exit(1);
    }
}

//======================================================================================
// CHECK WHICH OF THE ABOVE FUNCTIONS ARE PRIVATE!
private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceProperties& m_surf_props;
    WallState& m_wall_state;
    
    std::vector<WallProductionTerms*> v_surf_prod;

    DiffusionVelocityCalculator* mp_diff_vel_calc;
    MassBlowingRate* mp_mass_blowing_rate;

    // VARIABLES FOR SOLVER
    const size_t m_ns;
    double m_Twall;
    double m_Pwall;
    Eigen::VectorXd v_rhoi;
    Eigen::VectorXd v_work;
    Eigen::VectorXd v_X;
    Eigen::VectorXd v_dX;
    Eigen::VectorXd v_f;
    Eigen::MatrixXd v_jac;
    double m_pert;
    double m_X_unpert;
    Eigen::VectorXd v_f_unpert;

    const int set_state_with_rhoi_T;

    // TEMPORARY
    Eigen::VectorXd v_sep_mass_prod_rate;

};

//======================================================================================

Mutation::Utilities::Config::ObjectProvider<SurfaceBalanceSolverGamma, SurfaceBalanceSolver> surface_balance_solver_ablation("ablation"); // @totuesday
Mutation::Utilities::Config::ObjectProvider<SurfaceBalanceSolverGamma, SurfaceBalanceSolver> surface_balance_solver_gamma("gamma");
Mutation::Utilities::Config::ObjectProvider<SurfaceBalanceSolverGamma, SurfaceBalanceSolver> surface_balance_solver_frc("frc");
 
    } // namespace GasSurfaceInteraction
} // namespace Mutation
