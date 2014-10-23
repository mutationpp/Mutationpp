#include <cassert>

#include "GSIReaction.h"

using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

GSIReaction::GSIReaction(const IO::XmlElement& node, const Thermodynamics& thermo):
                      m_formula("")
{
    // Make sure this is a reaction type XML element
    assert( node.tag() == "reaction" );
}

void GSIReaction::parseFormula(const IO::XmlElement& node, const Thermodynamics& thermo)
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
    
//    if (m_formula[pos+1] == '>') {
//        m_reversible = false;
        products = m_formula.substr(pos+2, m_formula.length()-pos-1);
//    } else {
//        m_reversible = true;
//        products = m_formula.substr(pos+1, m_formula.length()-pos);
//    }
    

    // Now that we have reactant and product strings, we can parse each 
    // separately using the same algorithm
    parseSpecies(m_reactants, reactants, node, thermo);
    parseSpecies(m_products,  products,  node, thermo);
}

void GSIReaction::parseSpecies(
    std::vector<int>& species, std::string& str, const IO::XmlElement& node,
    const class Thermodynamics& thermo)
{  
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
            /*if (str.substr(s,e-s+1) == "M") {
                m_thirdbody = true;
            } else {*/
                int index = thermo.speciesIndex(str.substr(s,e-s+1));
                if (index >= 0)
                    for (int i = 0; i < nu; ++i)
                        species.push_back(index);
                else
                    node.parseError(("Species " + str.substr(s,e-s+1) +
                        " is not in the mixture list!").c_str());
            //}
            add_species = false;
        }
        
        // Move on to the next character
        c++;
    }
    
    // Sort the species vector
    std::sort(species.begin(), species.end()); 
}

/*
 *********************************************************************************************
 */

CatalysisReaction::CatalysisReaction(const IO::XmlElement& node, const Thermodynamics& thermo, const std::string m_category)
                       : GSIReaction(node, thermo)//,
                         //mp_catalysis_rate(NULL)
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
    IO::XmlElement::const_iterator iter = node.begin();
    for ( ; iter != node.end(); ++iter) {
        if (m_category == "gamma"){
            if (iter->tag() == "gamma_const"){
                //mp_catalysis_rate = new GammaModelConst(*iter);
            } else if(iter->tag() == "gamma_T"){
                /* Give an error here. Not implemented yet */;
            } else if(iter->tag() == "gamma_TP"){
                /* Give an error here. Not implemented yet */;
            } else {;/* Give an error here. Not implemented yet */}
        } else {/*Give an error here*/;}
    }   
}
