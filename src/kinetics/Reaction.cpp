#include <cstdlib>
#include <iostream>
#include <cassert>

#include "Reaction.h"

namespace Mutation {
    namespace Kinetics {

using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

void swap(Reaction& left, Reaction& right) {
    using std::swap;
    swap(left.m_formula,     right.m_formula);
    swap(left.m_reactants,   right.m_reactants);
    swap(left.m_products,    right.m_products);
    swap(left.m_reversible,  right.m_reversible);
    swap(left.m_thirdbody,   right.m_thirdbody);
    swap(left.m_thirdbodies, right.m_thirdbodies);
    swap(left.mp_rate,       right.mp_rate);
}

Reaction::Reaction(IO::XmlElement& node, const class Thermodynamics& thermo)
    : m_formula(""),
      m_reversible(true),
      m_thirdbody(false),
      mp_rate(NULL)
{
    // Make sure this is a reaction type XML element
    assert( node.tag() == "reaction" );

    // Store the reaction formula (must have)
    node.getAttribute("formula", m_formula, 
        "No formula specied with reaction!");
    
    // Parse the formula to determine which species are involved, whether or
    // not this is a third-body reaction, and reversibility of the reaction
    parseFormula(node, thermo);
    
    // Now loop through the children of this node to determine the other 
    // attributes of the reaction
    IO::XmlElement::Iterator iter = node.begin();
    for ( ; iter != node.end(); ++iter) {
        if (iter->tag() == "arrhenius") {
            mp_rate = new Arrhenius(*iter, order());
        } else if (iter->tag() == "M") {
            if (m_thirdbody) {
                std::vector<std::string> tokens;
                String::tokenize(iter->text(), tokens, ":, ");
                for (int i = 0; i < tokens.size(); i+=2) {
                    int index = thermo.speciesIndex(tokens[i]);
                    if (index >= 0) {
                        m_thirdbodies.push_back(
                            std::make_pair(index, atof(tokens[i+1].c_str())));
                    } else {
                        iter->parseError((
                            std::string("Thirdbody species ") + tokens[i] +
                            std::string(" is not in mixture list!")).c_str());
                    }
                }
            } else {
                iter->parseError((
                    std::string("This reaction is not a thirdbody reaction") +
                    std::string(" but thirdbodies are given!")).c_str());
            }
        }
    }
    
    // Make sure we got a RateLaw out of all that
    if (mp_rate == NULL)
        node.parseError("A rate law must be supplied with this reaction!");
    
    // Determine the reaction type
    determineType(thermo);
}

void determineType(const Thermodynamics& thermo)
{

}

void Reaction::parseFormula(
    IO::XmlElement& node, const class Thermodynamics& thermo)
{
    // First step is to split the formula into reactant and product
    // strings and determine reversibility of the reaction
    size_t pos = m_formula.find("=");
    if (pos == std::string::npos)
        node.parseError((
            std::string("Reaction formula ") + m_formula +
            std::string(" does not have '=' or '=>'!")).c_str());
    
    std::string reactants = m_formula.substr(0, pos);
    std::string products;
    
    if (m_formula[pos+1] == '>') {
        m_reversible = false;
        products = m_formula.substr(pos+2, m_formula.length()-pos-1);
    } else {
        m_reversible = true;
        products = m_formula.substr(pos+1, m_formula.length()-pos);
    }
    
    // Now that we have reactant and product strings, we can parse each 
    // separately using the same algorithm
    parseSpecies(m_reactants, reactants, node, thermo);
    parseSpecies(m_products,  products,  node, thermo);
}

void Reaction::parseSpecies(
    std::vector<int>& species, std::string& str, const IO::XmlElement& node,
    const class Thermodynamics& thermo)
{
    size_t c = 0;
    size_t s = 0;
    size_t e = 0;
    int nu   = 1;
    bool add_species = false;
    ParseState state = coefficient;
    
    // Remove all spaces to simplify the state machine logic
    String::eraseAll(str, " ");
    
    // Loop over every character
    while (c != str.size()) {
        // Decide what action to do
        switch(state) {
            case coefficient:
                if (isdigit(str[c])) {
                    nu = atoi(str.substr(c,1).c_str());
                    s = c + 1;
                } else {
                    nu = 1;
                    s = c;
                }                    
                state = name;
                break;
            case name:
                if (str[c] == '+')
                    state = plus;
                break;
            case plus:
                if (str[c] != '+') {
                    e = c - 2;                        
                    c--;
                    add_species = true;
                    state = coefficient;
                }
                break;                     
        }
        
        // Add the last species which would not be counted because c will
        // equal str.size() on the next iteration
        if (c == str.size() - 1) {
            add_species = true;
            e = c;
        }
        
        // If we found the start and end position of a species, add it to
        // the list nu times (unless it is the special case of 'M')
        if (add_species) {
            if (str.substr(s,e-s+1) == "M") {
                m_thirdbody = true;
            } else {
                int index = thermo.speciesIndex(str.substr(s,e-s+1));
                if (index >= 0)
                    for (int i = 0; i < nu; ++i)
                        species.push_back(index);
                else
                    node.parseError(("Species " + str.substr(s,e-s+1) +
                        " is not in the mixture list!").c_str());
            }
            add_species = false;
        }
        
        // Move on to the next character
        c++;
    }
}

    } // namespace Kinetics
} // namespace Mutation

