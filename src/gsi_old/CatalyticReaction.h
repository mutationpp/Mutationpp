#ifndef CATALYTIC_REACTION_H
#define CATALYTIC_REACTION_H

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
class CatalyticReaction{

public:
    CatalyticReaction( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, const Mutation::gsi::CatalysisSurfaceProperties* surf_props );
    ~CatalyticReaction();

    CatalyticReaction( const CatalyticReaction& l_catalytic_reaction )
               : m_formula( l_catalytic_reaction.m_formula ),
                 m_reactants( l_catalytic_reaction.m_reactants ),
                 m_products( l_catalytic_reaction.m_products ),
                 m_conserves( l_catalytic_reaction.m_conserves ),
                 m_has_active_sites( l_catalytic_reaction.m_has_active_sites ),
                 m_conserves_active_sites( l_catalytic_reaction.m_has_active_sites ),
                 mp_catalysis_rate( l_catalytic_reaction.mp_catalysis_rate != NULL ? l_catalytic_reaction.mp_catalysis_rate->clone() : NULL ){ }
    
    inline const std::string& formula() const { return m_formula; } 
    inline int nReactants() const { return m_reactants.size(); }
    inline int nProducts() const { return m_products.size(); }
    inline bool conservesChargeAndMass() const { return m_conserves; }
  
    const std::vector<int>& reactants() const { return m_reactants; }
    const std::vector<int>& products() const { return m_products;}

    void parseFormula( const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo );
    
    inline bool hasActiveSites() const { return m_has_active_sites; }
    inline bool conservesActiveSites() const { return m_conserves_active_sites; }
    const CatalysisRateLaw* CatalysisrateLaw() const { return mp_catalysis_rate; }

    void parseSpecies( std::vector<int>& species, std::string& str, const Mutation::Utilities::IO::XmlElement& node );
    
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;
    
    bool m_conserves;
    
private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const Mutation::gsi::CatalysisSurfaceProperties* mp_surf_props;
    CatalysisRateLaw* mp_catalysis_rate;
    bool m_has_active_sites;
    bool m_conserves_active_sites;
    
    void  errorInvalidRateLawGamma( const std::string& l_gamma_rate_law );
};

    } // namespace gsi
} //namespace Mutation

#endif // CATALYTIC_REACTION_H
