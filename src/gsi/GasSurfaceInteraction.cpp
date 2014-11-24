#include "GasSurfaceInteraction.h"
#include "Utilities.h"

/**
 * @todo list:
 * -> Make compilation of GSI obligatory. DONE
 * -> Add Ablation file. And everything that this includes.
 * -> Update the branch with trunk.
 * -> Finish with the SurfaceProperties parser.
 * -> Finish with the Rate Manager.
 * -> Finish with the Rate Laws.
 * -> Add the ability to understand "wall" phase species. Discuss about the generality of this implementation with the index.
 * -> In GSIReaction change thermo.nSpecies -> a member in order to reduce cost.
 */

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace gsi {
      
      using Mutation::Thermodynamics::Thermodynamics;
      
//==============================================================================
      
GasSurfaceInteraction::GasSurfaceInteraction(const Thermodynamics& thermo, std::string gsi_catalysis_mechanism_file/** @todo ablation , std::string gsi_ablation_mechanism_file*/):
                                              m_thermo(thermo)
{
    if (gsi_catalysis_mechanism_file == "none" /** @todo ablation && gsi_ablation_mechanism_file == "none"*/){ return;}
    
    if (gsi_catalysis_mechanism_file != "none"){
    
        gsi_catalysis_mechanism_file = 
        getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/gsi/" +
        gsi_catalysis_mechanism_file + ".xml";
    
        IO::XmlDocument doc(gsi_catalysis_mechanism_file);        
        IO::XmlElement root = doc.root();
    
        if (root.tag() != "gsi_mechanism") {
            std::cout << "Root element in gsi_mechanism file " << gsi_catalysis_mechanism_file
                      << " is not of 'gsi_mechanism' type!";
            exit(1); 
        }
        
        root.getAttribute("category", m_category, m_category);
        if (m_category != "catalysis") {
            std::cout << "The Gas-Surface interaction category in the " 
	              << gsi_catalysis_mechanism_file << " file is not of 'catalysis' type!";
            exit(1); 
        }
        
        /** @todo There are two ifs for the m_catalytic_model. One here and one in GSIReaction. Remove the other and only keep this one.*/
        root.getAttribute("model", m_catalytic_model, m_catalytic_model);
        if (m_catalytic_model == "gamma"){
	    /** @todo Here the rate manager should be added. */ 
            IO::XmlElement::const_iterator iter = root.begin();
            for ( ; iter != root.end(); ++iter) {        
                if (iter->tag() == "reaction")
                    addCatalyticReaction(CatalysisReaction(*iter, thermo, m_catalytic_model));
            }
        }
        else if (m_catalytic_model == "finite_rate_chemistry"){ 
	    /** @todo Here the rate manager should be added. */
            ;
        }
        else{
            std::cout << "The catalytic model " << m_catalytic_model 
            << " has not been implemented yet!";
            exit(1);  
        }
        
        // Setup the rate manager
        // mp_rates = new RateManager(thermo.nSpecies(), m_reactions);
        
    }       
    
    /** @todo ablation if (gsi_ablation_mechanism_file != "none"){ 
        std::cerr << "Ablation Module has not been implemented yet" << endl;
        exit(1);}*/
    
    //Closing the GSI reaction mechanisms
    closeGSIReactions(true);
        
}

//==============================================================================

GasSurfaceInteraction::~GasSurfaceInteraction()
{
    delete [] mp_rhoi;
    delete [] mp_Twall;
}


//==============================================================================

void GasSurfaceInteraction::netGSIProductionRates(double * const p_wdot)
{
  /**
   * @todo The general expression for the production rate for the gamma model it is 
   *
   *   omega_i = m_i M_i \sum_{reactions} \nu_ri * gamma_r - \sum_{r = reactions} \sum_{k = species} \mu_rik m_k M_k
   */
  
  
    // Special case of no reactions
//    if (nCatalyticReactions() == 0) {
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
             p_wdot[i] = 0.;
//    }
    
    
}

//==============================================================================
/**
 * @todo Passing a pointer. In the future maybe passing by reference, but create 
 *       the copy constructor.
 */

void GasSurfaceInteraction::addCatalyticReaction(const CatalysisReaction& catalytic_reaction)
{
    m_catalytic_reactions.push_back(catalytic_reaction);
}

//==============================================================================

//void GasSurfaceInteraction::addGSIReaction(const GSIReaction& gsi_reaction)
//{
//    m_gsi_reactions.push_back(gsi_reaction);
//}

//==============================================================================

void GasSurfaceInteraction::closeGSIReactions(const bool validate_gsi_mechanism)
{
   const size_t ns = m_thermo.nSpecies();
   
   // GSI Mechanism validation
   if (validate_gsi_mechanism) {
   
       bool is_valid = true;
       
       /**
	* @todo For this in the GSIReaction the Big Three should be implemented
	*  http://stackoverflow.com/questions/4172722/what-is-the-rule-of-three
	*/
       //Check for duplicate gsi reactions
       //for (size_t i = 0; i < nCatalyticReactions()-1; i++) {
       //    for (size_t j = i+1; j < nGSIReactions(); ++j) {
       //        if (m_gsi_reactions[i] == m_gsi_reactions)
       //    }
       //}
       
       // Check for elemental mass and charge conservation
//       for (size_t i = 0; i < nGSIReactions(); i++) {
//           if (!m_gsi_reactions[i].conservesChargeAndMass()){
//              std::cerr << "GSI reaction " << i+1 << " \"" << m_gsi_reactions[i].formula()
//                        << " \" does not conserve charge or mass." << endl;
//              is_valid = false;
//           }
//       }
       
//       // Check for the existance of active sites
//       for (size_t i = 0; i < nGSIReactions(); i++) {
//           if (!m_gsi_reactions[i].hasActiveSites()){
//              std::cerr << "GSI reaction " << i+1 << " \"" << m_gsi_reactions[i].formula()
//                        << " \" has no active sites or surface species." << endl;
//              is_valid = false;
//           }
//       }
       
       // Check for the conservation of active sites
//       for (size_t i = 0; i < nGSIReactions(); i++) {
//           if (!m_gsi_reactions[i].conservesActiveSites()){
//              std::cerr << "GSI reaction " << i+1 << " \"" << m_gsi_reactions[i].formula()
//                        << " \" does not conserve active sites" << endl;
//              is_valid = false;
//           }
//       }
       
       if (!is_valid) {
           cout << "Validation check of the gas surface interaction reactions failed!" << endl;
           exit(1);
       }
       
   }
   
}

//==============================================================================
      
void GasSurfaceInteraction::setWallState(const double* const p_rhoi, const double* const p_rhoie)
{
    /**
     * @todo Make it more general in case of multi temperature models.
     * @todo Overload it if needed
     * @todo Here I get a segmentation fault
     */
    mp_rhoi = new double (m_thermo.nSpecies());
    mp_Twall = new double (1);
    
    mp_Twall[0] = p_rhoie[0];
    
    
    for (int i = 0; i < m_thermo.nSpecies(); ++i) {
        mp_rhoi[i] = p_rhoi[i];
    }
    
}

      
    }// namespace GSI
} // namespace Mutation
