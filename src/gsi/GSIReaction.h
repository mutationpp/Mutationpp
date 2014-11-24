#ifndef GSI_REACTION_H
#define GSI_REACTION_H

#include "CatalysisRateLaws.h"
#include "Utilities.h"
#include "Thermodynamics.h"

namespace Mutation{
    namespace gsi{
/**
 * @todo Maybe make an abstract class GSIReaction with two children one Ablation and on Catalysis
**/
 
class GSIReaction
{
public:
    
    /**
     * Add description
     */
    GSIReaction(
        const Mutation::Utilities::IO::XmlElement& node, 
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    /**
     * Copy Constructor
     */
    GSIReaction(const GSIReaction& gsireaction)
        : m_formula(gsireaction.m_formula),
          m_reactants(gsireaction.m_reactants),
          m_products(gsireaction.m_products),
          m_conserves(gsireaction.m_conserves)
    { }
    
    
    /**
     * Destructor
     */
    virtual ~GSIReaction() { };
    
    /**
     * Returns the gsi reaction formula.
     */
    inline const std::string& formula() const {
        return m_formula;
    } 
    
    inline int nReactants() const {
        return m_reactants.size();
    }
    
    inline int nProducts() const {
        return m_products.size();
    }
  
    inline bool conservesChargeAndMass() const {
        return m_conserves;
    }
  
protected:
    /* virtual */ void parseFormula(const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo); // If it is needed
    
    /* virtual */ void parseSpecies(std::vector<int>& species,
        std::string& str,
        const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo); // If it is needed
    
protected:
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;
    
    bool m_conserves;
};

class CatalysisReaction: public GSIReaction{
public:

    CatalysisReaction(const Mutation::Utilities::IO::XmlElement& node,
		      const Mutation::Thermodynamics::Thermodynamics& thermo, 
		      const std::string m_category);
    
/** 
 * Copy constructor
 */
    
    CatalysisReaction(const CatalysisReaction& cat_reaction)
        : GSIReaction(cat_reaction),
          m_has_active_sites(cat_reaction.m_has_active_sites),
          m_conserves_active_sites(cat_reaction.m_has_active_sites),
          mp_catalysis_rate(cat_reaction.mp_catalysis_rate ? cat_reaction.mp_catalysis_rate->clone() : NULL)
    {}
    
    /**
     * Default Destructor
     */
    ~CatalysisReaction()
    {
        delete mp_catalysis_rate;
        mp_catalysis_rate = NULL;
    }
/**
 * 
 */
    inline bool hasActiveSites() const {
        return m_has_active_sites;
    }
  
/**
 * 
 */
    inline bool conservesActiveSites() const {
        return m_conserves_active_sites;
    }
    
    /**
     * 
     */
    const CatalysisRateLaw* CatalytisrateLaw() const { 
        return mp_catalysis_rate; 
    }
  
private:
    CatalysisRateLaw* mp_catalysis_rate;
    bool m_has_active_sites;
    bool m_conserves_active_sites;
    
};

// class AblationReaction: public GSIReaction{};

    } // namespace gsi
} //namespace Mutation

#endif // GSI_REACTION_H