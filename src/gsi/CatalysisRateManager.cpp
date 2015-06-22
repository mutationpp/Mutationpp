#include "CatalysisRateManager.h"

namespace Mutation{
    namespace gsi{
      
//=============================================================================
//=============================================================================
//=============================================================================

CatalysisGammaRateManager::CatalysisGammaRateManager( const Mutation::Thermodynamics::Thermodynamics& thermo,  const std::vector<CatalysisReaction>& catalytic_reactions)
        : m_ns( thermo.nSpecies() ),
          m_nr(catalytic_reactions.size()),
          m_reactions(catalytic_reactions),
          m_thermo(thermo)
{
    mp_gamma = new double [m_nr];
    mp_kf = new double [m_nr];

    for (int i_reac = 0; i_reac < m_nr; ++i_reac){
        m_reactants.addReaction( i_reac, m_reactions[i_reac].reactants());
        m_irr_products.addReaction(i_reac, m_reactions[i_reac].products());
    }
    //l_prod = new double [m_ns];
    mp_rop = new double [m_nr];

}
      
void CatalysisGammaRateManager::computeRate(const WallState& l_wall_state, double* const p_prod)
{
    for ( int i_reac = 0; i_reac < m_nr; ++i_reac ){
        mp_kf[i_reac] = m_reactions[i_reac].CatalysisrateLaw()->forwardReactionRate(l_wall_state.getWallRhoi(), l_wall_state.getWallTemperature()); 
    }

    std::fill(p_prod, p_prod+m_ns, 0.0);
    m_reactants.incrSpecies(mp_kf, p_prod);
    m_irr_products.decrSpecies(mp_kf, p_prod);

    for ( int i_species = 0; i_species < m_ns; ++i_species){
        p_prod[i_species] *= m_thermo.speciesMw(i_species);
    }

      /** @todo STRONG COMMENT ABOUT THE SIGN **/
}

//=============================================================================
//=============================================================================
//=============================================================================

CatalysisFRCRateManager::CatalysisFRCRateManager(const Mutation::Thermodynamics::Thermodynamics& thermo,
                             const Mutation::gsi::CatalysisSurfaceProperties* surf_props,
                             const std::vector<Mutation::gsi::CatalysisReaction>& reactions)
                             : m_ns(thermo.nSpecies()),
                               m_nr(reactions.size()),
                               m_nsites(surf_props->nSites()),
                               m_total_species_in_sites(surf_props->nTotalSpeciesinSites()),
                               m_reactions(reactions),
                               m_thermo(thermo)
{ 
    mp_species_in_site = new int [m_nsites];

    // Allocating memory for kf and kb;
    mp_kf = new double [m_nr];
    mp_kb = new double [m_nr]; // Basically I could make it that this memory is continuous

    mp_ropf = new double [m_nr]; // Make it better (NULL IT; DELETE IT CORRECTLY)
    mp_rop = new double [m_nr];

    for(int i_sites = 0; i_sites < m_nsites; ++i_sites ) {
        mp_species_in_site[i_sites] = surf_props->nSpeciesinSite(i_sites);
    }

    for(int i_reac = 0; i_reac < m_nr; i_reac++){
        m_reactants.addReaction(i_reac, m_reactions[i_reac].reactants());
        // if for rev
        m_irr_products.addReaction(i_reac, m_reactions[i_reac].products());
    }

   l_prod = new double [m_ns + m_total_species_in_sites + m_nsites](); // + m_n something else

}

//=============================================================================

void CatalysisFRCRateManager::computeRate(const WallState& l_wall_state, double* const p_prod) // Why not m_wall_state?
{

    // Computing species concentrations

    // Computing the forward reaction rates coefficients
      for ( int i_reac = 0; i_reac < m_nr; ++i_reac )
          mp_kf[i_reac] = m_reactions[i_reac].CatalysisrateLaw()->forwardReactionRate(l_wall_state.getWallRhoi(), l_wall_state.getWallTemperature()); //forwardReactionRate

    // Computing the backward reaction rate coefficient
//      for (int i_reac = 0; i_reac < nReactions(); ++i_reac )

    // Solving for steady state coverages at the wall
    solveSteadyStateCoverage(l_wall_state); // REMOVE THIS IN ORDER TO AVOID SEGMENTATION FAULT

    // Computing production rate
    for(int i_ns = 0; i_ns < m_ns; ++i_ns){
        p_prod[i_ns] = 0.E0;
    }

}

//=============================================================================

void CatalysisFRCRateManager::solveSteadyStateCoverage(const WallState& l_wall_state/** @todo is it needed?, double* const p_prod*/) // Something should be done here so that I do not loop a again over reactions and perform these calls to the function
{
  
    // Computing Functions
    //double* l_prod = new double [m_ns + m_total_species_in_sites]; // + m_n something else

    // Getting the initial configuration

    for(int i_species = 0 ; i_species < m_ns ; ++i_species){
       l_prod[i_species] = l_wall_state.getWallNumberDensities()[i_species];
    }

    for(int i_species_in_sites = 0 ; i_species_in_sites < m_total_species_in_sites + m_nsites ; ++i_species_in_sites)
        l_prod[i_species_in_sites + m_ns] = l_wall_state.getWallSiteDensities()[i_species_in_sites];

    for(int i_reac = 0; i_reac < m_nr; i_reac++)
        mp_ropf[i_reac] = mp_kf[i_reac];

    m_reactants.multReactions(l_prod, mp_ropf);

    for(int i_reac = 0; i_reac < m_nr; i_reac++)
        mp_rop[i_reac] = mp_ropf[i_reac];

    std::fill(l_prod, l_prod+(m_ns + m_total_species_in_sites + m_nsites), 0.0);
    m_reactants.decrSpecies(mp_rop, l_prod);
    m_irr_products.incrSpecies(mp_rop, l_prod);

    /** @todo IT SEEMS ALMOST CORRECT! CHECK IT AGAIN TOMORROW)

//    for(int i_ns = 0; i_ns < m_ns + m_total_species_in_sites + m_nsites ; i_ns++)
//        std::cout << "PROD = " << l_prod[i_ns] << std::endl;

    */

    // Computing Jacobians

    // Calling Newton Solver(Temporary A simple matrix inverter)

    // Saving Solution

    // Returning to computeRate through updating the initial solution

    // Updating the Surface Properties

    // delete [] l_prod; DELETE IT SOMEWHERE!!!

}


    } // namespace gsi
} // namespace Mutation
