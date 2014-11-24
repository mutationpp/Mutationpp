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
    GasSurfaceInteraction(const Mutation::Thermodynamics::Thermodynamics& thermo, const std::string gsi_catalysis_mechanism_file/** @todo ablation , std::string gsi_ablation_mechanism_file*/); 
    
/**
 * Destructor
 */
    ~GasSurfaceInteraction();
    
public:
  
/**
  * Returns the number of gsi reactions in the mechanism.
  */

    size_t nCatalyticReactions() const {
        m_catalytic_reactions.size();
    }

//    size_t nGSIReactions() const {
//        return m_gsi_reactions.size();
//    }
    
public:
/**
 * Add description
 */
    void netGSIProductionRates(double * const p_wdot);
   
public:
/**
 * Add description. Computes the Jacobian with respect to rho_i or y_i etc of (rho_i Vi) - wi = 0 
 */
    void jacobianRho();
    
public:
/**
 * Add description
 * @todo overload this function depending on what are the important conditions at the wall.
 */

    void setWallState(const double* const p_rhoi, const double* const p_rhoie);
    //void setWallState();
    
private:
/**
 * Add description
 */
    void addCatalyticReaction(const CatalysisReaction& catalytic_reaction);
    //void addAblativeReaction(const AblationReaction& ablative_reaction);
    //void addGSIReaction(const GSIReaction& gsi_reaction);
  
private:
/**
 * Add description
 */
    void closeGSIReactions(const bool validate_gsi_mechanism);

private:
/**
 * Private members
 */
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    
    //std::vector<GSIReaction> m_gsi_reactions;
    std::vector<CatalysisReaction> m_catalytic_reactions;
    //std::vector<AblationReaction> m_ablative_reactions;
    
    std::string m_category;
    std::string m_catalytic_model;
    
protected:
  
    /** 
     * @todo Fix this to be general, in accordance with StateModels
     */
    double* mp_Twall;
    double* mp_rhoi;
    
};

    } // namespace gsi
} //namespace Mutation

#endif // GSI_H