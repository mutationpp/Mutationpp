#ifndef GSI_H
#define GSI_H

#include <cassert>

#include "Utilities.h"
#include "Thermodynamics.h"
#include "CatalysisRateLaws.h"

namespace Mutation {
    namespace gsi {

/*
 ** General Variables
 */

//const Mutation::Thermodynamics::Thermodynamics* m_thermo;

/*
 ** Creates a structure with the important properties for the reactions. It will be a class soon. 
 ** Each reaction is saved in a vector.
 */

struct gsi_reaction{
  std::string m_formula;
  CatalysisRateLaw* mp_catalysis_rate;
  std::vector<int> m_reactants;
  std::vector<int> m_products;
  
  // Default Constructor
  gsi_reaction();
  
};

extern std::vector<gsi_reaction> m_gsi_reactions;

/*
 ** Initializer. First function to be called when initialing gsi 
 */

void gsi_initializer(const Mutation::Thermodynamics::Thermodynamics& thermo, std::string gsi_mechanism_file); 

/*
 ** Beginning with the parser
 */

gsi_reaction parser_wall_reactions(const Mutation::Utilities::IO::XmlElement& node); 
      
/*
 ** Function that adds reactions to a vector
 */
      
void add_gsi_reaction(const gsi_reaction& reaction);
      
/*
 ** Function that parses formula.
 */

void parse_formula(const Mutation::Utilities::IO::XmlElement& node, gsi_reaction& reaction); //const class Thermodynamics& thermo // gsi_reaction -> Probably removed when we have a class

/*
 ** Function that parses species
 */

//void parseSpecies(std::vector<int>& species, std::string& str, const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics* thermo); //const class Thermodynamics& thermo
void parseSpecies(std::vector<int>& species, std::string& str, const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics* thermo); 

    } // namespace GSI
} // namespace Mutation

#endif // GSI_H