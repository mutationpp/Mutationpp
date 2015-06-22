#include <cassert>

#include "GSIReaction.h"

using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

namespace Mutation{
    namespace gsi{
      
GSIReaction::GSIReaction(const IO::XmlElement& node, const class Thermodynamics& thermo, const CatalysisSurfaceProperties* surf_props):
                      m_formula(""),
                      m_conserves(true),
                      mp_surf_props(surf_props)
{
    // Make sure this is a reaction type XML element
    assert( node.tag() == "reaction" );
}

void GSIReaction::parseFormula(const IO::XmlElement& node, const class Thermodynamics& thermo)
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
    
    /** @todo Add reversible reactions **/
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
        int index;
        if (add_species) {
            index = thermo.speciesIndex(str.substr(s,e-s+1));
            
            if(index == -1){
                index = mp_surf_props->GSIspeciesIndex(str.substr(s,e-s+1));
            }
            
           if(index == -1){
                node.parseError(("Species " + str.substr(s,e-s+1) +
                    " is not in the mixture list or a species in the wall phase!").c_str());
            }

            species.push_back(index);
            add_species = false;
            
/*           if(e-s+1 >= 3){
                if(str.substr(e-1,2) == "-s")
                    index = thermo.speciesIndex(str.substr(s,e-s-1)) + thermo.nSpecies();
                //else if(str.substr(e-1,2) == "-b")
                //    index = thermo.speciesIndex(str.substr(s,e-s-1));
                else index = thermo.speciesIndex(str.substr(s,e-s+1));
            } else {
                if(str.substr(s,e-s+1) == "s") index = 2 * thermo.nSpecies();
                else index = thermo.speciesIndex(str.substr(s,e-s+1));
            }
            if (index >= 0)
                for (int i = 0; i < nu; ++i)
                    species.push_back(index);
            else
                node.parseError(("Species " + str.substr(s,e-s+1) +
                    " is not in the mixture list!").c_str());
            add_species = false; */

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

CatalysisReaction::CatalysisReaction(const IO::XmlElement& node, const class Thermodynamics& thermo, const std::string m_category, const CatalysisSurfaceProperties* surf_props)
                       : GSIReaction(node, thermo, surf_props), /** @todo Check again for the surf_props */
                         mp_catalysis_rate(NULL),
                         m_conserves_active_sites(true),
                         m_has_active_sites(false)
{
    // Make sure this is a reaction type XML element
    assert( node.tag() == "reaction" );
    
    // Store the reaction formula
    node.getAttribute("formula", m_formula, 
        "No formula specied with reaction!");
    
    // Parse the formula to determine which species are involved, whether or
    // not this is a third-body reaction, and reversibility of the reaction
    parseFormula(node, thermo);
  
    // Now loop through the children of this node to determine the other 
    // attributes of the reaction
    /**
     * @todo Here maybe no need to iterate
     */
    IO::XmlElement::const_iterator iter = node.begin();
    for ( ; iter != node.end(); ++iter) {
        if (m_category == "gamma"){
            if (iter->tag() == "gamma_const"){
                mp_catalysis_rate = new GammaModelConst(*iter, thermo, m_reactants);
            } else if(iter->tag() == "gamma_T"){
                std::cerr << "This category of gamma model, "           /** @todo */
                          << iter->tag() << ", has not been implemented yet!" << std::endl;      
                     exit(1);
            } else if(iter->tag() == "gamma_TP"){
                std::cerr << "This category of gamma model, "           /** @todo */
                          << iter->tag() << ", has not been implemented yet!" << std::endl;      
            } else {
                std::cerr << "This category of gamma model, "
                          << iter->tag() << ", has not been implemented yet!" << std::endl;      
            }
        } else if (m_category == "finite_rate_chemistry") {

            if (iter->tag() == "physisorption"){
                     mp_catalysis_rate = new Physisorption(*iter, thermo);
            } else if (iter->tag() == "thermal_desorption"){
                     mp_catalysis_rate = new ThermalDesorption(*iter, thermo);
            } else if(iter->tag() == "chemisorption") {
                     mp_catalysis_rate = new Chemisorption(*iter, thermo);
            } else if(iter->tag() == "E-R") {
                     mp_catalysis_rate = new ERRecombination(*iter,thermo);
            } else if(iter->tag() == "phys_to_chem") {
                     mp_catalysis_rate = new PhysisorptiontoChemisorption(*iter, thermo);
            } else if(iter->tag() == "L-H") {
                     mp_catalysis_rate = new LHRecombination(*iter, thermo);
            } else { std::cerr << "This mechanism of the Finite Rate Chemistry, " 
                     << iter->tag() << ", has not been implemented yet!" << std::endl;      
                     exit(1);
            }
        } else { std::cerr << "This category of Catalysis " << m_category << 
                 " has not been implemented yet!" << std::endl; 
                 exit(1);}
    }   
    
    // Make sure we got a RateLaw out of all that
    if (mp_catalysis_rate == NULL)
        node.parseError("A rate law must be supplied with this reaction!");
    
    // Check that reactions conserve elements and active sites
    /**
     * @todo To be rewritten
     */
//    const size_t ne = thermo.nElements();
//    int sums [ne];
//    std::fill(sums, sums+ne, 0);
//    int activesites = 0; // @todo: again here add the multiple active sites.
//    int forthermo;
//    for (int i = 0; i < nReactants(); ++i)
//        for (int k = 0; k < ne; ++k){
//            if (m_reactants[i] != 2 * thermo.nSpecies()){
//                if (m_reactants[i] >= thermo.nSpecies()){
//                    forthermo = m_reactants[i] - thermo.nSpecies();
//                    ++activesites;
//                    m_has_active_sites = true;
//                } else
//                    forthermo = m_reactants[i];
//                sums[k] += thermo.elementMatrix()(forthermo, k);
//            } else {
//                ++activesites;
//                m_has_active_sites = true;
//            }
//        }
//    for (int i = 0; i < nProducts(); ++i)
//        for (int k = 0; k < ne; ++k){
//            if (m_products[i] != 2 * thermo.nSpecies()){
//                if (m_products[i] >= thermo.nSpecies()){
//                    forthermo = m_products[i] - thermo.nSpecies();
//                   --activesites;
//                } else
//                    forthermo = m_products[i];
//                sums[k] -= thermo.elementMatrix()(forthermo, k);
//            } else {
//                --activesites;
//            }
//        }
//    for (int i = 0; i < ne; ++i){
//        m_conserves &= (sums[i] == 0);
//        m_conserves_active_sites &= (activesites == 0);
//    }
    
}

    } // namespace gsi
} //namespace Mutation
