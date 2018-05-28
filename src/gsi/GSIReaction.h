#ifndef GSI_REACTION_H
#define GSI_REACTION_H

#include "DataGSIReaction.h"

namespace Mutation {
    namespace GasSurfaceInteraction {
/**
 * Abstract base class that defines a complete reaction. Contains information
 * for the reactants, products, reversibility and rate laws. Solid classes can
 *  be of type catalysis, ablation and FRC. Species information is stored via
 *  indices into the global species list and surface properties species.
 */
class GSIReaction
{
public:
	/**
	 * Required for self registering different types of reactions.
	 */
    typedef const DataGSIReaction& ARGS;

    /*
     * Constructor of the abstract base class, taking as input a structure
     * of type DataGSIReaction
     */
    GSIReaction(ARGS args)
        : mp_rate_law(NULL){ }

	/// Returns name of this type.
	static std::string typeName() { return "GSIReaction"; }

    /**
     * Destructor
     */
    virtual ~GSIReaction(){}

    /**
     * Returns pointer to the rate law object associated with this reaction.
     */
    GSIRateLaw* getRateLaw() const { return mp_rate_law; }

    /**
     * Returns a constant reference to the vector of the reactants in this
     * reactions
     */
    const std::vector<int>& getReactants() const { return m_reactants; }

    /**
     * Returns a constant reference to the vector of the products in this
     * reactions
     */
    const std::vector<int>& getProducts() const { return m_products; }

protected:
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;

    GSIRateLaw* mp_rate_law;

    /**
     * Error function if no reaction formula has been provided.
     */
    inline const char* errorNoFormulainReaction() const {
    	return "No formula specied with reaction!"; }

    /**
     * Determines the reactants and products of a heterogeneous reaction
     * given the reaction formula.  The formula should be in the form
     * "aA + bB => cC + dD" where (a,b,c,d) are the optional stoichiometric
     * coefficients (A,B,C,D) represent species names. Spaces and number of
     * species are completely optional. The formula can have up to three
     * reactants and three products.
     */
    void parseFormula(
        const Mutation::Thermodynamics::Thermodynamics& l_thermo,
        const Mutation::Utilities::IO::XmlElement& l_node_reaction,
        const SurfaceProperties& l_surf_props)
    {
        std::string l_reactants;
        std::string l_products;
        splitFormulainReactantsProducts(l_reactants, l_products, l_node_reaction);
        
        parseSpecies(m_reactants, l_reactants, l_node_reaction, l_thermo, l_surf_props);
        parseSpecies(m_products, l_products, l_node_reaction, l_thermo, l_surf_props);
    }
    
    /**
     * Splits reaction formula string from the XmlElement in a string of
     * reactants and products in order for them to be parsed into species.
     */
    virtual void splitFormulainReactantsProducts(
        std::string& l_reactants, std::string& l_products,
        const Mutation::Utilities::IO::XmlElement& l_node_reaction)
    {
        size_t l_pos_equal = m_formula.find("=");
        if(l_pos_equal == std::string::npos){
            l_node_reaction.parseError(
                (std::string("Reaction formula ") + m_formula
                 + std::string(" does not have '=' or '=>'!")).c_str());
        }

        l_reactants = m_formula.substr(0, l_pos_equal);
        l_products = m_formula.substr(l_pos_equal + 2, m_formula.length() - l_pos_equal - 1);
    }

    /**
     * Purely virtual class which fills a list of strings with the species
     * names of the species in a reaction formula with the format
     * "aA + bB + cC" where (a,b,c) represents the stoichiometric coefficients
     * of the (A,B,C) species names.  The species can be present either in the
     * gas or in the surface phase. The coefficients must be positive integers
     * in the range 1-9 and species names can have a '+' at the end of the
     * name so this is taken into account.  In addition, the coefficients
     * default to 1 if they are not present in front of the species name. For
     * coefficients other than 1, the species is added multiple times to the
     * species list.
     */
    virtual void parseSpecies(
        std::vector<int>& l_species, std::string& l_str_chem_species,
        const Mutation::Utilities::IO::XmlElement& l_node_reaction,
        const Mutation::Thermodynamics::Thermodynamics& l_thermo,
        const SurfaceProperties& l_surf_props) = 0;
}; // class GSIReaction

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSI_REACTION_H
