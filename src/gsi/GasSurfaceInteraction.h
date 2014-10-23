#ifndef GSI_H
#define GSI_H

#include "Utilities.h"
#include "Thermodynamics.h"
#include "GSIReaction.h"


namespace Mutation {
    namespace gsi {
      
class GasSurfaceInteraction
{
  
public:
/**
 * Add description
 */
    GasSurfaceInteraction(const Mutation::Thermodynamics::Thermodynamics& thermo, std::string gsi_catalysis_mechanism_file/** @todo ablation , std::string gsi_ablation_mechanism_file*/); 

private:
  /**
   * Add description
   */
    void addGSIReaction(const GSIReaction &gsi_reaction);
  
  
private:
  /**
   * Private members
   */
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    
    std::vector<GSIReaction> m_gsi_reactions;
    
    std::string m_category;
    std::string m_catalytic_model;
  
  
  
};

    } // namespace gsi
} //namespace Mutation

#endif // GSI_H