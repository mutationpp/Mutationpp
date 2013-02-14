#ifndef KINETICS_REACTION_H
#define KINETICS_REACTION_H

#include <string>
#include <vector>
#include <set>

#include "RateLaws.h"

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
    Reaction(Mutation::Utilities::IO::XmlElement& node);
    
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
    int reactant(const std::string& species) const {
        return m_reactants.count(species);
    }
    
    /**
     * Returns the multiset of species acting as reactants in this reaction.
     */
    const std::multiset<std::string>& reactants() const {
        return m_reactants;
    }
    
    /**
     * Returns the product stoichiometric coefficient for the given species.
     */
    int product(const std::string& species) const {
        return m_products.count(species);
    }
    
    /**
     * Returns the multiset of species acting as products in this reaction.
     */
    const std::multiset<std::string>& products() const {
        return m_products;
    }
    
    /**
     * Returns a vector of (species name, thirdbody efficiency) value pairs.  If
     * this reaction is not a thirdbody reaction then the vector will be empty.
     * Note that only species with non-default efficiencies are listed.
     */
    const std::vector<std::pair<std::string, double> >& efficiencies() const {
        return m_thirdbodies;
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
    void parseFormula(Mutation::Utilities::IO::XmlElement& node);
    
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
    void parseSpecies(std::multiset<std::string>& species, std::string& str);

private:

    std::string m_formula;
    std::multiset<std::string> m_reactants;
    std::multiset<std::string> m_products;
    
    bool m_reversible;
    bool m_thirdbody;
    
    std::vector<std::pair<std::string, double> > m_thirdbodies;
    
    RateLaw* mp_rate;

}; // class Reaction


/**
 * Swaps one reaction for the other and vice-versa.
 */
void swap(Reaction& left, Reaction& right);


#endif // KINETICS_REACTION_H
