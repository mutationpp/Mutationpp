#include "GasSurfaceInteraction.h"
#include "Utilities.h"

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace gsi {
      
      using Mutation::Thermodynamics::Thermodynamics;
      
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
        
        root.getAttribute("model", m_catalytic_model, m_catalytic_model);
        if (m_catalytic_model == "gamma"){
	  
            IO::XmlElement::const_iterator iter = root.begin();
            for ( ; iter != root.end(); ++iter) {        
            if (iter->tag() == "reaction")
                addGSIReaction(CatalysisReaction(*iter, thermo, m_catalytic_model));
            }
            
        }
            
        }
        else{
            std::cout << "The catalytic model " << m_catalytic_model 
            << " has not been implemented yet!";
            exit(1);  
  
        }
        

        
    /** @todo ablation if (gsi_ablation_mechanism_file != "none"){ 
        std::cerr << "Ablation Module has not been implemented yet" << endl;
        exit(1);}*/
        
}

void GasSurfaceInteraction::addGSIReaction(const GSIReaction& gsi_reaction)
{
m_gsi_reactions.push_back(gsi_reaction);
}


      
    }// namespace GSI
} // namespace Mutation
