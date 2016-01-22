#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIReaction.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReactionFRC : public GSIReaction {

public:
    GSIReactionFRC( ARGS l_data_gsi_reaction )
                     : GSIReaction( l_data_gsi_reaction ) {

        assert( l_data_gsi_reaction.s_iter_reaction->tag() == "reaction" );

        l_data_gsi_reaction.s_iter_reaction->getAttribute( "formula", m_formula, errorNoFormulainReaction() );
        parseFormula( l_data_gsi_reaction.s_thermo, *(l_data_gsi_reaction.s_iter_reaction), l_data_gsi_reaction.s_surf_props );

        const Mutation::Utilities::IO::XmlElement& l_node_rate_law = *(l_data_gsi_reaction.s_iter_reaction->begin());
        DataGSIRateLaw l_data_gsi_rate_law = { l_data_gsi_reaction.s_thermo,
                                               l_node_rate_law, 
                                               l_data_gsi_reaction.s_surf_props,
                                               m_reactants };
        mp_rate_law = Mutation::Utilities::Config::Factory<GSIRateLaw>::create( l_node_rate_law.tag(), l_data_gsi_rate_law );

        if ( mp_rate_law == NULL ) {
            l_data_gsi_reaction.s_iter_reaction->parseError( "A rate law must be provided for this reaction!" );
        }

    }

//=============================================================================================================

    ~GSIReactionFRC(){ 
        if ( mp_rate_law != NULL ){ delete mp_rate_law; }
    }

//=============================================================================================================

    bool isCatalytic(){ return 0; }
    bool isAblative(){ return 0; }

//=============================================================================================================

    void parseSpecies( std::vector<int>& l_species, std::string& l_str_chem_species, 
                       const Mutation::Utilities::IO::XmlElement& l_node_reaction, 
                       const Mutation::Thermodynamics::Thermodynamics& l_thermo,
                       const SurfaceProperties& l_surf_props ){
        // State-Machine states for parsing a species formula
        enum ParseState {
            coefficient,
            name,
            plus
        } state = coefficient;

        size_t c = 0;
        size_t s = 0;
        size_t e = 0;
        int nu   = 1;
        bool add_species = false;
        
        // Remove all spaces to simplify the state machine logic
        Mutation::Utilities::String::eraseAll(l_str_chem_species, " ");
        
        // Loop over every character
        while (c != l_str_chem_species.size()) {
            // Decide what action to do
            switch(state) {
                case coefficient:
                    if (isdigit(l_str_chem_species[c])) {
                        nu = atoi(l_str_chem_species.substr(c,1).c_str());
                        s = c + 1;
                    } else {
                        nu = 1;
                        s = c;
                    }                    
                    state = name;
                    break;
                case name:
                    if (l_str_chem_species[c] == '+')
                        state = plus;
                    break;
                case plus:
                    if (l_str_chem_species[c] != '+') {
                        e = c - 2;                        
                        c--;
                        add_species = true;
                        state = coefficient;
                    }
                    break;                     
            }
            
            // Add the last species which would not be counted because c will
            // equal str.size() on the next iteration
            if (c == l_str_chem_species.size() - 1) {
                add_species = true;
                e = c;
            }
            
            // If we found the start and end position of a species, add it to
            // the list nu times (unless it is the special case of 'M')
            int index;
            if (add_species) {
                index = l_thermo.speciesIndex(l_str_chem_species.substr(s,e-s+1));
                
                if(index == -1){
                    index = l_surf_props.speciesIndexWall(l_str_chem_species.substr(s,e-s+1));
                }

               if(index == -1){
                    l_node_reaction.parseError(("Species " + l_str_chem_species.substr(s,e-s+1) +
                        " is not in the mixture list or a species in the wall phase!").c_str());
                }

                l_species.push_back(index);
                add_species = false;
            }
            
            // Move on to the next character
            c++;
        }
        
        // Sort the species vector
        std::sort(l_species.begin(), l_species.end()); 
        }

//=============================================================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIReactionFRC, GSIReaction> frc_reaction("frc");

    } // namespace GasSurfaceInteraction
}  // namespace Mutation
