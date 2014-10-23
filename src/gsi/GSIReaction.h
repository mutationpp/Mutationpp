#ifndef GSI_REACTION_H
#define GSI_REACTION_H

#include "Utilities.h"
#include "Thermodynamics.h"
#include "CatalysisRateLaws.h"

//Maybe make an abstract class GSIReaction with two children one Ablation and on Catalysis

class GSIReaction
{
public:
    
    /**
     * Add description
     */
    GSIReaction(
        const Mutation::Utilities::IO::XmlElement& node, 
        const Mutation::Thermodynamics::Thermodynamics& thermo);
  
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
};

class CatalysisReaction: public GSIReaction{
public:

    CatalysisReaction(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, const std::string m_category);
  
private:
    //CatalysisRateLaw* mp_catalysis_rate;
    
};

// class AblationReaction: public GSIReaction{};

#endif // GSI_REACTION_H