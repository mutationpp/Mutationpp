#include "CatalysisRateManager.h"

namespace Mutation{
    namespace gsi{
      
//=============================================================================
//=============================================================================
//=============================================================================

CatalysisRateManagerGamma::CatalysisRateManagerGamma( const Mutation::Thermodynamics::Thermodynamics& thermo,  const std::vector<CatalysisReaction>& catalytic_reactions)
        : m_ns( thermo.nSpecies() ),
          v_wall_reaction_rate_constant( thermo.nSpecies() ),
          m_nr(catalytic_reactions.size()),
          v_reactions(catalytic_reactions),
          m_thermo(thermo){

    for (int i_reac = 0; i_reac < m_nr; ++i_reac){
        v_reactants.addReaction( i_reac, v_reactions[i_reac].reactants() );
        m_irr_products.addReaction( i_reac, v_reactions[i_reac].products() );
    }

}
      
//=============================================================================

void CatalysisRateManagerGamma::computeRate(const WallState& l_wall_state, double* const lp_mass_production_rate){
 
    // Get reaction rate constant
    for ( int i_reac = 0; i_reac < m_nr; ++i_reac ){
        v_wall_reaction_rate_constant(i_reac) = v_reactions[i_reac].CatalysisrateLaw()->forwardReactionRate(l_wall_state.getWallRhoi(), l_wall_state.getWallTemperature()); 
    }

    // Constant rate times 
    std::fill(lp_mass_production_rate, lp_mass_production_rate + m_ns, 0.0);
    v_reactants.incrSpecies( &v_wall_reaction_rate_constant(0), lp_mass_production_rate );
    m_irr_products.decrSpecies( &v_wall_reaction_rate_constant(0), lp_mass_production_rate );

    // Multiply by molar mass
    for ( int i_species = 0; i_species < m_ns; ++i_species){
        lp_mass_production_rate[i_species] *= m_thermo.speciesMw(i_species);
    }

}

//=============================================================================
//=============================================================================
//=============================================================================

CatalysisRateManagerFiniteRateChemistry::CatalysisRateManagerFiniteRateChemistry(const Mutation::Thermodynamics::Thermodynamics& thermo,
                             const Mutation::gsi::CatalysisSurfaceProperties* surf_props,
                             const std::vector<Mutation::gsi::CatalysisReaction>& lv_reactions)
                             : m_ns( thermo.nSpecies() ),
                               m_nr ( lv_reactions.size() ),
                               m_nsites( surf_props->nSites()),
                               m_total_species_in_sites( surf_props->nTotalSpeciesinSites() ),
                               v_reactions( lv_reactions ),
                               v_wall_reaction_rate_constant( lv_reactions.size() ),
                               v_species_in_site( surf_props->nSites() ),
                               v_mass_production_rate( thermo.nSpecies() + surf_props->nTotalSpeciesinSites() + surf_props->nSites() ),
                               m_thermo(thermo){
 
    

    mp_ropf = new double [m_nr]; 
    mp_rop = new double [m_nr];

    for( int i_sites = 0; i_sites < m_nsites; ++i_sites ) {
        v_species_in_site[i_sites] = surf_props->nSpeciesinSite(i_sites);
    }

    for( int i_reac = 0; i_reac < m_nr; i_reac++){
        v_reactants.addReaction(i_reac, v_reactions[i_reac].reactants());
        m_irr_products.addReaction(i_reac, v_reactions[i_reac].products());
    }

}

//=============================================================================

CatalysisRateManagerFiniteRateChemistry::~CatalysisRateManagerFiniteRateChemistry(){ 

        delete [] mp_ropf;
        delete [] mp_rop;
    }

//=============================================================================

void CatalysisRateManagerFiniteRateChemistry::computeRate(const WallState& l_wall_state, double* const lp_mass_production_rate) // Why not m_wall_state?
{

    // Computing species concentrations

    // Computing the forward reaction rates coefficients
    for ( int i_reac = 0; i_reac < m_nr; ++i_reac ){
        v_wall_reaction_rate_constant(i_reac) = v_reactions[i_reac].CatalysisrateLaw()->forwardReactionRate( l_wall_state.getWallRhoi(), l_wall_state.getWallTemperature()); //forwardReactionRate
    }

    // Solving for steady state coverages at the wall
    solveSteadyStateCoverage( l_wall_state ); 

    // Computing production rate
    for(int i_ns = 0; i_ns < m_ns; ++i_ns){
        lp_mass_production_rate[i_ns] = 0.E0;
    }

}

//=============================================================================

void CatalysisRateManagerFiniteRateChemistry::solveSteadyStateCoverage( const WallState& l_wall_state ) 
{
  
//        // Computing Functions
//        double* l_prod = new double [m_ns + m_total_species_in_sites]; // + m_n something else
//    
//        // Getting the initial configuration
//    
//    //    l_prod->v_mass_production_rate;
//    
//        for(int i_species = 0 ; i_species < m_ns ; ++i_species){
//           l_prod[i_species] = l_wall_state.getWallNumberDensities()[i_species];
//        }
//    
//        for(int i_species_in_sites = 0 ; i_species_in_sites < m_total_species_in_sites + m_nsites ; ++i_species_in_sites)
//            l_prod[i_species_in_sites + m_ns] = l_wall_state.getWallSiteDensities()[i_species_in_sites];
//    
//        for(int i_reac = 0; i_reac < m_nr; i_reac++)
//            mp_ropf[i_reac] = v_wall_reaction_rate_constant(i_reac);
//    
//        v_reactants.multReactions(l_prod, mp_ropf);
//    
//        for(int i_reac = 0; i_reac < m_nr; i_reac++)
//            mp_rop[i_reac] = mp_ropf[i_reac];
//    
//        std::fill(l_prod, l_prod+(m_ns + m_total_species_in_sites + m_nsites), 0.0);
//        v_reactants.decrSpecies(mp_rop, l_prod);
//        m_irr_products.incrSpecies(mp_rop, l_prod);
//    
//        // Computing Jacobians
//    
//        // Saving Solution
//    
//        // Returning to computeRate through updating the initial solution
//    
//        // Updating the Surface Properties
//    
//        delete [] l_prod; 

}


    } // namespace gsi
} // namespace Mutation
