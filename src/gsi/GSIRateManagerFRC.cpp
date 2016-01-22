#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"

#include <Eigen/Dense>

#include "NewtonSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerFRC : public GSIRateManager, public Mutation::Numerics::NewtonSolver<Eigen::VectorXd, GSIRateManagerFRC> {

public:
    GSIRateManagerFRC( DataGSIRateManager l_data_rate_manager ) 
                     : GSIRateManager( l_data_rate_manager ) {

    	for ( int i_reac = 0 ; i_reac < l_data_rate_manager.s_reactions.size() ; ++i_reac ){
    	    m_reactants.addReaction(i_reac, l_data_rate_manager.s_reactions[i_reac]->getReactants() );
            m_irr_products.addReaction( i_reac, l_data_rate_manager.s_reactions[i_reac]->getProducts() );
    	}

    }

//=============================================================================

    ~GSIRateManagerFRC(){ }

//=============================================================================

    void computeRate( Eigen::VectorXd& lv_mass_prod_rate ){

    	Eigen::VectorXd v_full_prod_rate( m_thermo.nSpecies() + m_surf_props.nSpeciesWall() );
    	// ! v_conc_wall(  m_thermo.nSpecies() + m_surf_props.nSpeciesWall() )
    	// Setting up the surface state in the full vector
    	m_wall_state.getConcWallSurfState( v_full_prod_rate );

    	// Getting the kfs with the initial conditions
    	Eigen::VectorXd v_wall_reac_rate( v_reactions.size() ); // Make this global in the class!
        for ( int i_reac = 0; i_reac < v_wall_reac_rate.size() ; ++i_reac ) {
            v_wall_reac_rate(i_reac) =  v_reactions[i_reac]->
            		getRateLaw()->forwardReactionRate( m_wall_state.getWallRhoi(),
            				                           m_wall_state.getWallT()    );
        }

        // Solving for the wall species
//        computeSurfSpeciesBalance( v_initial_wall_state );

        // Update the wall species in wall state

        // Getting the gas phase production rates in mol
        m_reactants.incrSpecies( v_wall_reac_rate, v_full_prod_rate );
        m_irr_products.decrSpecies( v_wall_reac_rate, v_full_prod_rate );

        // Multiply by molar mass
        for ( int i_sp = 0 ; i_sp < lv_mass_prod_rate.size() ; ++i_sp ){
            lv_mass_prod_rate(i_sp) = v_full_prod_rate(i_sp) * m_thermo.speciesMw(i_sp);
        }

    }

//=============================================================================

private:

void computeSurfSpeciesBalance(){

    // solve()

	// setWallStateSurf

}
//=============================================================================

//    Eigen::VectorXd v_wall_reac_rate;

//=============================================================================

// Wall Steady State Solver
    void updateFunction( Eigen::VectorXd& v_X ){ // Takes as input the current state of the wall!

    	// updateStateVector // of wall and of gas
    	// v_F.block(0, 0, 0, n_sp ) = v_mol_wall;
    	// for ns->ns+n_surf
    	// v_F.block(0, n_sp ) = v_X;

    	// Initialize v_fun to 0.0
    	// m_reactants.incrSpeciesSurface( v_wall_reac_rate, v_tota_X );
    	// m_irr_products.decrSpeciesSurface( v_wall_reac_rate, v_tota_X );

    	// v_fun is updated
    	for ( int i_surf_species = 0 ; i_surf_species < m_surf_props.nSpeciesWall() ; ++i_surf_species ){
    	    v_F( i_surf_species ) =  ( i_surf_species + m_thermo.nSpecies() ) ;
    	}

    }

    void updateJacobian( Eigen::VectorXd& v_X ){

    	// Save Unperturbed solution
    	v_F_unper = v_F;


        // Make pert
        //    updateFunction();
    	// Fill in correct v_jac places


    	// Unpert
    	v_F = v_F_unper;

    }

    double norm(){ v_F.lpNorm<Eigen::Infinity>(); }
    Eigen::VectorXd systemSolution(){ return v_jac.partialPivLu().solve( v_F ); }

    Eigen::VectorXd v_F; // v_F
    Eigen::VectorXd v_F_unper; // v_F_unper
    Eigen::MatrixXd v_jac; // v_J

    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateManagerFRC, GSIRateManager> gsi_rate_manager_frc("frc");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
