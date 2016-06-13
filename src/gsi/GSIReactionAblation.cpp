#include "GSIReaction.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReactionAblation : public GSIReaction
{
public:
    GSIReactionAblation(ARGS args)
        : GSIReaction(args)
    {
        assert(args.s_iter_reaction->tag() == "reaction");

        args.s_iter_reaction->getAttribute("formula", m_formula, errorNoFormulainReaction());
        parseFormula(args.s_thermo, *(args.s_iter_reaction), args.s_surf_props);

        const Mutation::Utilities::IO::XmlElement& l_node_rate_law = *(args.s_iter_reaction->begin());
        DataGSIRateLaw l_data_gsi_rate_law = { args.s_thermo,
                                               args.s_transport,
                                               l_node_rate_law,
                                               args.s_surf_props,
                                               m_reactants,
                                               m_products };

        mp_rate_law = Mutation::Utilities::Config::Factory<GSIRateLaw>::create(
            l_node_rate_law.tag(), l_data_gsi_rate_law);

        if (mp_rate_law == NULL) {
            args.s_iter_reaction->parseError("A rate law must be provided for this reaction!");
        }
    }

//=============================================================================================================

    ~GSIReactionAblation(){
        if (mp_rate_law != NULL){ delete mp_rate_law; }
    }

//=============================================================================================================

    void parseSpecies(
        std::vector<int>& l_species, std::string& l_str_chem_species,
        const Mutation::Utilities::IO::XmlElement& l_node_reaction,
        const Mutation::Thermodynamics::Thermodynamics& l_thermo,
        const SurfaceProperties& l_surf_props){

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
        
        Mutation::Utilities::String::eraseAll(l_str_chem_species, " ");
        
        // Loop over every character
        while (c != l_str_chem_species.size()) {
            switch(state) {
                case coefficient:
                    if (isdigit(l_str_chem_species[c])) {
                        nu = atoi(l_str_chem_species.substr(c, 1).c_str());
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
            if (c == l_str_chem_species.size() - 1) {
                add_species = true;
                e = c;
            }

            int index;
            if (add_species) {
                index = l_thermo.speciesIndex(l_str_chem_species.substr(s, e-s+1));
                // Check if the parsed input species are solid carbon.
                // To be FIXED. The surface species should exist in the surface
                // properties!

                if(l_str_chem_species.substr( s, e-s+1 ) == "C-s"){ // @toGeorge
                	index = -2;
                }

                if(index == -1){
                    l_node_reaction.parseError(( "Species " + l_str_chem_species.substr( s, e-s+1 )
                                                  + " is not in the mixture list or a species in the wall phase!" ).c_str() );
                }

                if (index != -2){
                    l_species.push_back(index);
                }

                add_species = false;
            }

            // Move on the next character
            c++;
        }

        // Sort the species vector
        std::sort(l_species.begin(), l_species.end());
    }
}; //class GSIReactionAblation

Mutation::Utilities::Config::ObjectProvider<GSIReactionAblation, GSIReaction> ablation_reaction("ablation");

    } // namespace GasSurfaceInteraction
}  // namespace Mutation
