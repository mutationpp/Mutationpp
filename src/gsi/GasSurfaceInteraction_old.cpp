
#include "GasSurfaceInteraction.h"

using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

namespace Mutation {
    namespace gsi {
      
using Mutation::Thermodynamics::Thermodynamics;

/*
 ***************************************************************************************************
*/
const Mutation::Thermodynamics::Thermodynamics* m_thermo;

/*
 ***************************************************************************************************
*/

gsi_reaction::gsi_reaction()
    : m_formula(""),
      mp_catalysis_rate(NULL)
{}
  
/*
 ***************************************************************************************************
*/

std::vector<gsi_reaction> m_gsi_reactions;
    
/*
 ***************************************************************************************************
*/

//void gsi_initializer(Thermodynamics& thermo, std::string gsi_mechanism)
void gsi_initializer(const Thermodynamics& thermo, std::string gsi_mechanism)
{
    m_thermo = &thermo;
  
    std::string m_category;
    
    if (gsi_mechanism == "none")
    return;
    
    gsi_mechanism = 
    getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/gsi/" +
    gsi_mechanism + ".xml";
    
    IO::XmlDocument doc(gsi_mechanism);        
    IO::XmlElement root = doc.root();
    
    if (root.tag() != "gsi_mechanism") {
        std::cout << "Root element in gsi_mechanism file " << gsi_mechanism
             << " is not of 'gsi_mechanism' type!";
        exit(1); 
    }
    
    root.getAttribute("category", m_category, m_category);
    if (m_category == "catalysis"){
    // Now loop over all of the reaction nodes and add each reaction to the
    // corresponding data structure pieces
    IO::XmlElement::const_iterator iter = root.begin();
    for ( ; iter != root.end(); ++iter) {        
        if (iter->tag() == "reaction")
            add_gsi_reaction(parser_wall_reactions(*iter));
    }
    }
    else if(m_category == "ablation"){;
    }
    else{;}
    
  
}

/*
 ***************************************************************************************************
*/

gsi_reaction parser_wall_reactions(const IO::XmlElement& node){
  
    gsi_reaction local_struc_gsi_reaction;
  
    // Make sure this is a reaction type XML element
    assert( node.tag() == "reaction" );
    
    // Store the reaction formula (must have)
    node.getAttribute("formula", local_struc_gsi_reaction.m_formula, 
        "No formula specied with reaction!");
    
    // Parse the formula to determine which species are involved, whether or
    // not this is a third-body reaction, and reversibility of the reaction
    parse_formula(node, local_struc_gsi_reaction);
    
    // HERE WE SHOULD HAVE A CHECK WHETHER IT IS, CATALYSIS ABLATION OR MIX
    
    // Now loop through the children of this node to determine the other 
    // attributes of the reaction
    IO::XmlElement::const_iterator iter = node.begin();
    for ( ; iter != node.end(); ++iter) {
        if (iter->tag() == "gamma_const"){
            local_struc_gsi_reaction.mp_catalysis_rate = new GammaModelConst(*iter);
        } else if(iter->tag() == "gamma_T"){
	    ;
	} else if(iter->tag() == "gamma_TP"){
	    ;
        } else {;}
    }   
    
    // Make sure we got a RateLaw out of all that
    if (local_struc_gsi_reaction.mp_catalysis_rate == NULL)
    node.parseError("A catalytic model must be supplied with this reaction!");
    
    // Check for charge and mass conservation
/*    const size_t ne = thermo.nElements();
    int sums [ne];
    std::fill(sums, sums+ne, 0.0);
    for (int i = 0; i < nReactants(); ++i)
        for (int k = 0; k < ne; ++k)
            sums[k] += thermo.elementMatrix()(m_reactants[i],k);
    for (int i = 0; i < nProducts(); ++i)
        for (int k = 0; k < ne; ++k)
            sums[k] -= thermo.elementMatrix()(m_products[i],k);
    for (int i = 0; i < ne; ++i)
        m_conserves &= (sums[i] == 0);
    
    // Figure out what type of reaction this is
    //determineType(thermo);*/
    return(local_struc_gsi_reaction);
  
}

/*
 ***************************************************************************************************
*/

void add_gsi_reaction(const gsi_reaction& reaction)
{
    m_gsi_reactions.push_back(reaction);
}


/*
 ***************************************************************************************************
*/

void parse_formula(const IO::XmlElement& node, gsi_reaction& local_struc_gsi_reaction)
{
    // First step is to split the formula into reactant and product
    // strings and determine reversibility of the reaction
    size_t pos = local_struc_gsi_reaction.m_formula.find("=");
    if (pos == std::string::npos)
        node.parseError((
            std::string("Reaction formula ") + local_struc_gsi_reaction.m_formula +
            std::string(" does not have '=' or '=>'!")).c_str());
    
    std::string reactants = local_struc_gsi_reaction.m_formula.substr(0, pos);
    std::string products;
    
//    if (m_formula[pos+1] == '>') {
//        m_reversible = false;
        products = local_struc_gsi_reaction.m_formula.substr(pos+2, local_struc_gsi_reaction.m_formula.length()-pos-1);
//    } else {
//        m_reversible = true;
//        products = m_formula.substr(pos+1, m_formula.length()-pos);
//    }
    

    // Now that we have reactant and product strings, we can parse each 
    // separately using the same algorithm
    parseSpecies(local_struc_gsi_reaction.m_reactants, reactants, node, m_thermo);
    parseSpecies(local_struc_gsi_reaction.m_products,  products,  node, m_thermo);

}

/*
 ***************************************************************************************************
*/

//void parseSpecies(std::vector< int >& species, std::string& str, const IO::XmlElement& node, Mutation::Thermodynamics::Thermodynamics* thermo)
void parseSpecies(std::vector<int>& species, std::string& str, const IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics* thermo) 
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
                int index = thermo->speciesIndex(str.substr(s,e-s+1));
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
 ***************************************************************************************************
*/
    } // namespace gsi
} // namespace Mutation
