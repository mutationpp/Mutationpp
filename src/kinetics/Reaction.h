#ifndef KINETICS_REACTION_H
#define KINETICS_REACTION_H

#include <string>
#include <vector>

#include "Thermodynamics.h"
#include "RateLaws.h"

using Mutation::Thermodynamics::ELECTRON;

namespace Mutation {
    namespace Kinetics {
    
enum ReactionType
{
    DISSOCIATION,
    EXCHANGE,
    IONIZATION,
    RECOMBINATION,
    ELECTRONIC_ATTACHMENT,
    ASSOCIATIVE_IONIZATION,
    CHARGE_EXCHANGE,
    ION_RECOMBINATION,
    ELECTRONIC_DETACHMENT,
    DISSOCIATIVE_RECOMBINATION
};

/**
 * Stores information that defines a complete reaction (reactants, products,
 * reversibility, thirdbody efficiencies, rate law, and rate coefficients).  The
 * Reaction has no knowledge of the available species in a mixture and species
 * information is stored via the species names, not actual Species objects.
 */
class Reaction
{
public:

    /**
     * Constructs a reaction object from a "reaction" XML element.
     */
    Reaction(
        Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Copy constructor.
     */
    Reaction(const Reaction& reaction)
        : m_formula(reaction.m_formula),
          m_reactants(reaction.m_reactants),
          m_products(reaction.m_products),
          m_reversible(reaction.m_reversible),
          m_thirdbody(reaction.m_thirdbody),
          m_thirdbodies(reaction.m_thirdbodies),
          m_type(reaction.m_type),
          mp_rate(reaction.mp_rate ? reaction.mp_rate->clone() : NULL)
    { }
    
    /**
     * Assignment operator.
     */
    Reaction& operator = (Reaction reaction) {
        swap(*this, reaction);
        return *this;
    }
    
    /**
     * Destructor.
     */
    ~Reaction() {
        delete mp_rate;
        mp_rate = NULL;
    }
    
    /**
     * Returns the reaction formula.
     */
    const std::string& formula() const {
        return m_formula;
    }
    

    /**
    * Returns the type of reaction.
    */
    const ReactionType type() const {
        return m_type;
    }
    
    /**
     * Returns the order of the reaction.  Note that this can be different from
     * the number of reactants in the reaction due to thirdbodies.
     */
    int order() const { 
        return m_reactants.size() + (m_thirdbody ? 1 : 0); 
    }
    
    /**
     * Returns the number of reactants (not including thirdbodies).
     */
    int nReactants() const { 
        return m_reactants.size(); 
    }
    
    /**
     * Returns the number of products (not including thirdbodies).
     */
    int nProducts() const { 
        return m_products.size(); 
    }
    
    /**
     * Returns true if the reaction is reversible, false otherwise.
     */
    bool isReversible() const { 
        return m_reversible; 
    }
    
    /**
     * Returns true if the reaction is a thirdbody reaction, false otherwise.
     */
    bool isThirdbody() const { 
        return m_thirdbody; 
    }
    
    /**
     * Returns the rate law object associated with this reaction.
     */
    const RateLaw* rateLaw() const { 
        return mp_rate; 
    }
    
    /**
     * Returns the reactant stoichiometric coefficient for the given species.
     */
    int reactant(int species) const {
        int count = 0;
        for (int i = 0; i < nReactants(); ++i)
            count += (m_reactants[i] == species ? 1 : 0);
        return count;
    }
    
    /**
     * Returns the multiset of species acting as reactants in this reaction.
     */
    const std::vector<int>& reactants() const {
        return m_reactants;
    }
    
    /**
     * Returns the product stoichiometric coefficient for the given species.
     */
    int product(int species) const {
        int count = 0;
        for (int i = 0; i < nReactants(); ++i)
            count += (m_products[i] == species ? 1 : 0);
        return count;
    }
    
    /**
     * Returns the multiset of species acting as products in this reaction.
     */
    const std::vector<int>& products() const {
        return m_products;
    }
    
    /**
     * Returns a vector of (species name, thirdbody efficiency) value pairs.  If
     * this reaction is not a thirdbody reaction then the vector will be empty.
     * Note that only species with non-default efficiencies are listed.
     */
    const std::vector<std::pair<int, double> >& efficiencies() const {
        return m_thirdbodies;
    }
    
    /**
    * Returns true when ions are present in the reactants
    **/
    bool ionReactants(const Mutation::Thermodynamics::Thermodynamics& thermo) const {
        for (int i=0; i < nReactants(); ++i) 
        if (thermo.species(m_reactants[i]).isIon() && thermo.species(m_reactants[i]).type() != ELECTRON) return true;
        return false;
    }

    /**
    * Returns true when ions are present in the products
    **/
    bool ionProducts(const Mutation::Thermodynamics::Thermodynamics& thermo) const {
        for (int i=0; i < nProducts(); ++i) 
           if (thermo.species(m_products[i]).isIon()) return true;
        return false;
    }

    /**
    * Returns true when electrons in reactants
    **/
    bool electronReactants(const Mutation::Thermodynamics::Thermodynamics& thermo) const {
        for (int i=0; i < nReactants(); ++i) 
           if (thermo.species(m_reactants[i]).type()==ELECTRON) return true;
        return false;
    }

    /**
    * Returns true when electrons in products
    **/
    bool electronProducts(const Mutation::Thermodynamics::Thermodynamics& thermo) const {
        for (int i=0; i < nProducts(); ++i) 
           if (thermo.species(m_products[i]).type()==ELECTRON) return true;
        return false;
    }


    friend void swap(Reaction& left, Reaction& right);
    
private:

    /**
     * Determines the reactants, products, reversibility, and whether or not
     * the reaction is a thirdbody reaction given the reaction formula.  The 
     * formula should be in the form "aA + bB + M = cC + dD" where (a,b,c,d) are
     * the optional stoichiometric coefficients (A,B,C,D) represent species
     * names.  Spaces and number of species are completely optional but there
     * should be at least one reactant and one product.  The '=' denotes a 
     * reversible reaction while '=>' could be used to denote a irreversible
     * forward reaction.
     */
    void parseFormula(
        Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Enumerates the states used in the state machine which parses the
     * reactants and products from the reaction formula.
     * @see parseSpecies()
     */
    enum ParseState {
        coefficient,
        name,
        plus
    };
    
    /**
     * Fills a list of strings with the species names of the species in a
     * reaction formula with the format "aA + bB + cC" where (a,b,c) represents
     * the stoichiometric coefficients of the (A,B,C) species names.  The 
     * coefficients must be positive integers in the range 1-9 and species names 
     * can have a '+' at the end of the name so this is taken into account.  In 
     * addition, the coefficients default to 1 if they are not present in front 
     * of the species name.  For coefficients other than 1, the species is added 
     * multiple times to the species list. 
     */
    void parseSpecies(
        std::vector<int>& species, std::string& str,
        const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);

    /**
    *Determines reaction type.  
    */
    void determineType(
       const Mutation:: Thermodynamics::Thermodynamics& thermo);

private:

    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;
    
    bool m_reversible;
    bool m_thirdbody;
    
    std::vector<std::pair<int, double> > m_thirdbodies;
    
    ReactionType m_type;
    RateLaw* mp_rate;

}; // class Reaction


/**
 * Swaps one reaction for the other and vice-versa.
 */
void swap(Reaction& left, Reaction& right);


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_REACTION_H
