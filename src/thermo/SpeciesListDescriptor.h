#ifndef THERMO_SPECIES_DESCRIPTOR_H
#define THERMO_SPECIES_DESCRIPTOR_H

#include <string>
#include <list>
#include <vector>
#include <set>

namespace Mutation {
    namespace Thermodynamics {

class Species;

/**
 * This class is used by thermodynamic databases to decide which species they 
 * should load upon initialization.  The species list could be determined from
 * a simple list of species names, or something more complex such as all gases
 * containing certain elements.
 */
class SpeciesListDescriptor
{
public:
    
    /**
     * Constructor which takes a string representation of the descriptor.
     */
    SpeciesListDescriptor(std::string descriptor);
    
    /**
     * Destructor.
     */
    ~SpeciesListDescriptor() { }
    
    /**
     * Tests if a species object is described by this descriptor and keeps track
     * of whether or not all explicitly defined species are matched.
     */
    bool matches(const Species& species) const;
    
    /**
     * Orders the species given as input in the output array.  Ensures that 
     * species explicitly listed by the user maintain the same order and that 
     * the electron is always at the beginning if it is present as a species.
     * Condensed phase species are listed at the end of the array.
     */
    void order(
        std::list<Species>& input, std::vector<Species>& output,
        std::vector<std::string>& missing) const;
    
private:

    /// Explicitly defined species names
    std::vector<std::string> m_species_names;
    
    /// List of allowed elements
    std::set<std::string> m_element_names;
    
    /// List of species who should be expanded into excited electronic states
    std::set<std::string> m_expand_states;
    
    /// True if gases are allowed
    bool m_gases;
    
    /// True if solids are allowed
    bool m_solids;
    
    /// True if liquids are allowed
    bool m_liquids;
    
};

    } // namespace Thermodynamics
} // namespace Mutation

#endif

