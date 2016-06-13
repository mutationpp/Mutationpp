#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"

#include <Eigen/Dense>

#include "NewtonSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerFRC : public GSIRateManager, public Mutation::Numerics::NewtonSolver<Eigen::VectorXd, GSIRateManagerFRC>
{
public:

    GSIRateManagerFRC(DataGSIRateManager args)
        : GSIRateManager(args),
          m_ns_gas(args.s_thermo.nSpecies()),
          m_ns_wall(args.s_surf_props.nSpeciesWall()),
          m_nr(v_reactions.size()),
          v_kf(m_nr),
          v_RateofProduction(m_nr),
          lv_numb_dens(m_ns_gas + m_ns_wall),
          v_X(m_ns_wall),
          v_F(m_ns_wall),
          v_F_unper(m_ns_wall),
          v_jac(m_ns_wall, m_ns_wall),
          v_F_total(m_ns_gas+m_ns_wall),
          m_pert(1.e-4) // @todo tolerance!
    {
    	for (int i_reac = 0 ; i_reac < m_nr ; ++i_reac){
    	    m_reactants.addReaction(i_reac, args.s_reactions[i_reac]->getReactants());
            m_irr_products.addReaction(i_reac, args.s_reactions[i_reac]->getProducts());
    	}
    }

//=============================================================================

    ~GSIRateManagerFRC(){}

//=============================================================================

    Eigen::VectorXd computeRate()
    {
    	// Getting the initial number densities
    	m_wall_state.getConcWallSurfState(lv_numb_dens); // @todo return the value

    	// Getting the kfs with the initial conditions
        for (int i_reac = 0; i_reac < m_nr; ++i_reac) {
            v_kf(i_reac) =  v_reactions[i_reac]->
            		            getRateLaw()->forwardReactionRateCoefficient(
            		                m_wall_state.getWallRhoi(), m_wall_state.getWallT());
        }

        // Solving for the steady state for the Surface Species!
        computeSurfSpeciesBalance();

        v_RateofProduction = v_kf;
        m_reactants.multReactions(lv_numb_dens, v_RateofProduction);

        lv_numb_dens.setZero();
        m_reactants.decrSpecies(v_RateofProduction, lv_numb_dens);
        m_irr_products.incrSpecies(v_RateofProduction, lv_numb_dens);

        // Multiply by molar mass
        return (lv_numb_dens.cwiseProduct(
        		Eigen::Map<const Eigen::VectorXd>(m_thermo.speciesMw(), m_ns_gas))/Mutation::NA).head(m_ns_gas);
    }

//=============================================================================

private:
    void computeSurfSpeciesBalance()
    {
      // Getting Number Densities!
      v_X.head( m_surf_props.nSpeciesWall() ) = lv_numb_dens.tail( m_surf_props.nSpeciesWall() );
      v_X = solve(v_X);
      lv_numb_dens.tail( m_surf_props.nSpeciesWall() ) = v_X.head( m_surf_props.nSpeciesWall() );

/*   @todo Make sure to keep l_total_n_sites
     @todo setWallStateSurf
*/
    }
//=============================================================================

public:
// Wall Steady State Solver
    void updateFunction(Eigen::VectorXd& lv_X)
    {
    	v_F_total.head(m_ns_gas) = lv_numb_dens.head(m_ns_gas);
    	v_F_total.tail(m_ns_wall) = lv_X;

        v_RateofProduction = v_kf;
        m_reactants.multReactions(v_F_total, v_RateofProduction);

        v_F_total.setZero();
    	m_reactants.decrSpecies(v_RateofProduction, v_F_total);
    	m_irr_products.incrSpecies(v_RateofProduction, v_F_total);

    	v_F = v_F_total.tail(m_ns_wall) ;

    	// Getting the conservation of Sites // THERE (IS)WAS A BUG
    	int pos_in_v_F = 0;
        int pos_in_v_X = 0;
    	double n_total_sites = m_surf_props.nTotalSites();
        int n_sp_in_site;
        double l_frac_site;
        for (int i_surf = 0; i_surf < m_surf_props.nSites() ; ++i_surf){
    		n_sp_in_site = m_surf_props.nSpeciesSite(i_surf);
    		l_frac_site = m_surf_props.fracSite(i_surf);
    		v_F(pos_in_v_F) = n_total_sites * l_frac_site;
            for(int i_sp_site = 0; i_sp_site < n_sp_in_site; i_sp_site++) {
                v_F(pos_in_v_F) -= lv_X(pos_in_v_X);
                pos_in_v_X++;
            }
            pos_in_v_F += n_sp_in_site;
    	}
    }

    void updateJacobian(Eigen::VectorXd& lv_X){

    	// Save Unperturbed solution
    	v_F_unper = v_F;

        // Make pert
    	for (int i = 0; i < m_ns_wall; ++i){
    		m_X_unper = lv_X(i);
    		lv_X(i) += lv_X(i) * m_pert;
            updateFunction( lv_X );
            v_jac.col(i) = (v_F - v_F_unper)/ (m_pert * m_X_unper) ;
            lv_X(i) = m_X_unper;
    	}

    	// Unpert
    	v_F = v_F_unper;

    }

    double norm(){ return v_F.lpNorm<Eigen::Infinity>(); }
    Eigen::VectorXd systemSolution(){ return v_jac.partialPivLu().solve(v_F); }

private:

    const size_t m_ns_gas;
    const size_t m_ns_wall;
    const size_t m_nr;

    // Solver Properties
    Eigen::VectorXd v_kf;
    Eigen::VectorXd v_RateofProduction;

    Eigen::VectorXd v_X; // v_X
    Eigen::VectorXd v_F; // v_F
    Eigen::VectorXd v_F_unper; // v_F_unper
    Eigen::MatrixXd v_jac; // v_J

    Eigen::VectorXd lv_numb_dens;
    Eigen::VectorXd v_F_total;

    double m_pert;
    double m_X_unper;

    // GSIStoichiometryManager
    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;

}; // class GSIRateManagerFRC

Mutation::Utilities::Config::ObjectProvider<GSIRateManagerFRC, GSIRateManager> gsi_rate_manager_frc("frc");

    } // namespace GasSurfaceInteraction
} // namespace Mutation


/*  Debugging Input...
 *
     	lv_numb_dens(2) /= 2;
    	lv_numb_dens(3) *= 2;
    	lv_numb_dens(4) /= 2;
    	lv_numb_dens(5) *= 2;
    	std::cout << "Number Densities= " << std::endl;
        std::cout << lv_numb_dens << std::endl;
 *
        std::cout << "Printing kfs = " << std::endl;
        std::cout << v_kf << std::endl;

        //  The following lines are for DEBUGGING purposes ONLY. REMOVE THEM
        m_reactants.multReactions( lv_numb_dens, v_RateofProduction );  // numb_dens, kf
        std::cout << "AFTER BEGIN" << std::endl;
    	std::cout << "Number Densities= " << std::endl;
        std::cout << lv_numb_dens << std::endl;
        std::cout << "Printing kfs = " << std::endl;
        std::cout << v_RateofProduction << std::endl;
        std::cout << "AFTER END" << std::endl; // TESTED TILL HERE

        lv_numb_dens.setZero();
        m_reactants.decrSpecies( v_RateofProduction, lv_numb_dens );
        std::cout << "AFTER AFTER BEGIN" << std::endl;
    	std::cout << "Number Densities= " << std::endl;
        std::cout << lv_numb_dens << std::endl;
        std::cout << "Printing kfs = " << std::endl;
        std::cout << v_RateofProduction << std::endl;
        std::cout << "AFTER AFTER END" << std::endl;

        m_irr_products.incrSpecies( v_RateofProduction, lv_numb_dens );
        std::cout << "AFTER AFTER AFTER BEGIN" << std::endl;
    	std::cout << "Number Densities= " << std::endl;
        std::cout << lv_numb_dens << std::endl;
        std::cout << "Printing kfs = " << std::endl;
        std::cout << v_RateofProduction << std::endl;
        std::cout << "AFTER AFTER AFTER END" << std::endl;
        std::cout << "Production Rates in #/m3 = " << std::endl;
        std::cout << lv_numb_dens << std::endl;
        // TILL HERE
*/
