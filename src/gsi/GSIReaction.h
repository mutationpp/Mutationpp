#ifndef GSIREACTION_H
#define GSIREACTION_H

#include <string>
#include <vector>

#include "DataGSIReaction.h"
#include "GSIRateLaw.h"
#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction {

public:
    typedef const DataGSIReaction& ARGS;

    GSIReaction( ARGS l_data_gsi_reaction )
               : mp_rate_law( NULL ){ }
    virtual ~GSIReaction(){ }

    GSIRateLaw* getRateLaw() const { return mp_rate_law; }

    const std::vector<int>& getReactants() const { return m_reactants; }
    const std::vector<int>& getProducts() const { return m_products; }

protected:
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;

    GSIRateLaw* mp_rate_law;

    inline const char* errorNoFormulainReaction() const { return "No formula specied with reaction!"; }

    void parseFormula( const Mutation::Thermodynamics::Thermodynamics& l_thermo, 
                       const Mutation::Utilities::IO::XmlElement& l_node_reaction,
                       const SurfaceProperties& l_surf_props ){

        std::string l_reactants;
        std::string l_products;
        splitFormulainReactantsProducts( l_reactants, l_products, l_node_reaction );
        
        parseSpecies( m_reactants, l_reactants, l_node_reaction, l_thermo, l_surf_props );
        parseSpecies( m_products, l_products, l_node_reaction, l_thermo, l_surf_props );

    }
    
    virtual void splitFormulainReactantsProducts( std::string& l_reactants, std::string& l_products, 
                                                  const Mutation::Utilities::IO::XmlElement& l_node_reaction ){

        size_t l_pos_equal = m_formula.find( "=" );
        if( l_pos_equal == std::string::npos ){
            l_node_reaction.parseError( ( std::string( "Reaction formula " ) + m_formula
                                        + std::string( " does not have '=' or '=>'!" ) ).c_str() );
        } // @todo test

        l_reactants = m_formula.substr( 0, l_pos_equal );
        l_products = m_formula.substr( l_pos_equal + 2, m_formula.length() - l_pos_equal - 1 );

    }

    virtual void parseSpecies ( std::vector<int>& l_species, std::string& l_str_chem_species, 
                                const Mutation::Utilities::IO::XmlElement& l_node_reaction, 
                                const Mutation::Thermodynamics::Thermodynamics& l_thermo,
                                const SurfaceProperties& l_surf_props ) = 0;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSIREACTION_H
