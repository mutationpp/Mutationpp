/**
 * @file Reaction.h
 *
 * @brief Declaration of Reaction class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef KINETICS_REACTION_H
#define KINETICS_REACTION_H

#include <string>
#include <vector>

#include "Thermodynamics.h"
#include "RateLaws.h"
#include "ReactionType.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Stores information that defines a complete reaction (reactants, products,
 * reversibility, thirdbody efficiencies, rate law, and rate coefficients).
 * Species information is stored via indices into the global species list.
 */
class Reaction
{
public:

    /**
     * Constructs a reaction object from a "reaction" XML element.
     */
    Reaction(
        const Mutation::Utilities::IO::XmlElement& node,
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
          m_conserves(reaction.m_conserves),
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
     * Equal to operator.  Returns true if these reactions have the same
     * stoichiometry in either the forward or reverse directions.
     */
    bool operator == (const Reaction& r);
    
    /**
     * Not equal to operator.
     */
    bool operator != (const Reaction& r) {
        return !((*this) == r);
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
     * Returns true if the reaction conserves mass and charge.
     */
    bool conservesChargeAndMass() const {
        return m_conserves;
    }
    
    /**
     * Returns the type of the reaction.
     * @see enum ReactionType.
     */
     ReactionType type() const {
        return m_type;
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
        const Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Enumerates the states used in the state machine which parses the
     * reactants and products from the reaction formula.
     * @see parseSpecies()
     */
    /*enum ParseState {
        coefficient,
        name,
        plus
    };*/
    
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
     * Determines which type of reaction this is.
     * @see enum ReactionType.
     */
    void determineType(const Mutation::Thermodynamics::Thermodynamics& thermo);

private:

    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;
    
    bool m_reversible;
    bool m_thirdbody;
    bool m_conserves;
    
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
