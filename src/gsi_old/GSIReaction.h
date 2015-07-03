#ifndef GSI_REACTION_H
#define GSI_REACTION_H

#include <string>

#include "CatalysisRateLaws.h"
#include "SurfaceProperties.h"
#include "Utilities.h"
#include "Thermodynamics.h"

namespace Mutation{
    namespace gsi{

/**
 * Base class from which the specific reactions for catalysis and ablation are derived.
 */
class GSIReaction
{
public:
    GSIReaction( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, const Mutation::gsi::CatalysisSurfaceProperties* surf_props );

    GSIReaction( const GSIReaction& gsireaction )
               : m_formula( gsireaction.m_formula ),
                 m_reactants( gsireaction.m_reactants ),
                 m_products( gsireaction.m_products ),
                 m_conserves( gsireaction.m_conserves ){ }
    
    virtual ~GSIReaction() { }
    
    inline const std::string& formula() const { return m_formula; } 
    inline int nReactants() const { return m_reactants.size(); }
    inline int nProducts() const { return m_products.size(); }
    inline bool conservesChargeAndMass() const { return m_conserves; }
  
    const std::vector<int>& reactants() const { return m_reactants; }
    const std::vector<int>& products() const { return m_products;}

protected:
    void parseFormula( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo );
    
    void parseSpecies(std::vector<int>& species,
    std::string& str,
    const Mutation::Utilities::IO::XmlElement& node,
    const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;
    
    bool m_conserves;
    
private:
    const Mutation::gsi::CatalysisSurfaceProperties* mp_surf_props;
};

class CatalysisReaction: public GSIReaction{
public:

    CatalysisReaction( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, const std::string m_category, const Mutation::gsi::CatalysisSurfaceProperties* surf_props );
    
    CatalysisReaction( const CatalysisReaction& cat_reaction )
                     : GSIReaction( cat_reaction ),
                       m_has_active_sites( cat_reaction.m_has_active_sites ),
                       m_conserves_active_sites( cat_reaction.m_has_active_sites ),
                       mp_catalysis_rate( cat_reaction.mp_catalysis_rate != NULL ? cat_reaction.mp_catalysis_rate->clone() : NULL ){ }
    
    ~CatalysisReaction();
 
    inline bool hasActiveSites() const { return m_has_active_sites; }
    inline bool conservesActiveSites() const { return m_conserves_active_sites; }
    const CatalysisRateLaw* CatalysisrateLaw() const { return mp_catalysis_rate; }
  
private:
    CatalysisRateLaw* mp_catalysis_rate;
    bool m_has_active_sites;
    bool m_conserves_active_sites;
    
    void  errorInvalidRateLawGamma( const std::string& l_gamma_rate_law );

};

    } // namespace gsi
} //namespace Mutation

#endif // GSI_REACTION_H
