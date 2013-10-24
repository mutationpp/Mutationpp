#ifndef THERMO_SPECIES_H
#define THERMO_SPECIES_H

#include <map>
#include <string>
#include <set>
#include <vector>

// Forward declaration of the XmlElement Class
namespace Mutation {
namespace Utilities {
namespace IO {
    class XmlElement;
} } }
     

namespace Mutation {
    namespace Thermodynamics {
    
class ParticleRRHO; 

/**
 * Responsible for loading element data from the elements.xml file.
 */
class Element
{
public:
    
    /**
     * Loads an element from an "element" XmlElement node.
     */
    Element(Mutation::Utilities::IO::XmlElement& xml_element);
    
    /**
     * Returns the name of the element.
     */
    const std::string &name() const {
        return m_name;
    }
    
    /**
     * Returns the atomic mass of the element in kg/mol.
     */
    double atomicMass() const { 
        return m_atomic_mass; 
    }
    
    /**
     * Returns the charge of the element.
     */
    int charge() const { 
        return m_charge; 
    } 
    
private:

    std::string m_name;
    double m_atomic_mass;
    int m_charge;

}; // class Element

/**
 * Enumerates the different phase types that are possible.
 */
enum PhaseType {
    GAS = 0,
    LIQUID,
    SOLID
};

/**
 * Enumerates the different particle types that are possible.
 */
enum ParticleType {
    ELECTRON,
    ATOM,
    MOLECULE
};

/**
 * Responsible for loading species data from the species.xml file.
 */
class Species
{
public:
    
    /**
     * Load species data from a "species" XmlElement node.
     */
    Species(
        Mutation::Utilities::IO::XmlElement &xml_element, 
        const std::vector<Element> &elements, std::set<int> &used_elements);
    
    /**
     * Instantiate a new species object which represents a single electronic
     * state of the given species.
     */
    Species(const Species& species, const size_t level);
    
    /**
     * Copy constructor.
     */
    Species(const Species& to_copy);
    
    /**
     * Destructor.
     */
    ~Species();
    
    /**
     * Assignment operator.
     */
    Species& operator = (Species species) {
        swap(*this, species);
        return *this;
    }
    
    /**
     * Returns the species name.
     */
    const std::string& name() const { 
        return m_name; 
    }
    
    /**
     * Returns the name of this species that is used in the NASA polynomial
     * database given the number of coefficients of the polynomial (7 or 9).
     */
    const std::string& nasaName(int n) const 
    {
        switch (n) {
            case 7:
                return m_nasa7_name;
            case 9:
                return m_nasa9_name;
            default:
                return m_name;
        }
    }
    
    /**
     * Returns type of species (ie: electron, atom, or molecule).
     *
     * @see enum ParticleType
     */
    ParticleType type() const {
        return (name() == "e-" ? ELECTRON : (nAtoms() == 1 ? ATOM : MOLECULE));
    }
    
    /**
     * Returns the molecular weight of the species in kg/mol.
     */
    const double &molecularWeight() const { return m_mw; }
    
    /**
     * Returns the number of atoms of a particular element contained in this
     * species.  For example calling nAtoms("C") for the C2H2 molecule would
     * return 2.
     */
    const int &nAtoms(const std::string &element) {
        return m_stoichiometry[element];
    }
    
    /**
     * Returns the total number of atoms in this species.  For example, C2H2 has
     * 4 total atoms (2 Carbon and 2 Hydrogen atoms).  Note that this does not
     * include electrons in the calculation.
     */
    int nAtoms() const {
        int atoms = 0;
        std::map<std::string, int>::const_iterator iter = 
            m_stoichiometry.begin();
        
        for ( ; iter != m_stoichiometry.end(); ++iter)
            atoms += (iter->first != "e-" ? iter->second : 0);
        
        return atoms;
    }
    
    /**
     * Returns the charge of this species.
     */
    int charge() const { 
        return m_charge; 
    }
    
    /**
     * Returns true if this is an ion.
     */
    bool isIon() const {
        return (charge() != 0);
    }
    
    /**
     * Returns the phase of this species.
     */
    PhaseType phase() const { 
        return m_phase; 
    }
    
    /**
     * Returns true if this species data contains a Rigid-Rotator Harmonic-
     * Oscillator model.
     */
    bool hasRRHOParameters() const {
        return (mp_rrho_model != NULL);
    }
    
    /**
     * Returns a pointer to this species ParticleRRHO object which stores
     * parameters for evaluating the species thermodynamics using the Rigid-
     * Rotator Harmonic-Oscillator model.  Note that if hasRRHOParameters() 
     * returns false, then this method will return a NULL pointer.
     */
    const ParticleRRHO* getRRHOParameters() const {
        return mp_rrho_model;
    }
    
    friend void swap(Species&, Species&);

private:

    /**
     * Determines the stoichiometry of this species from a stoichiometry string.
     */
    void loadStoichiometry(
        const std::string& stoichiometry, const std::vector<Element> &elements);

    /**
     * Will check the species name against its stoichiometry to determine if
     * the stoichiometry given by the user matches the stoichiometry found by
     * parsing the name.  Note that the name can be anything, therefor if it is
     * not in the standard form then it will simply be ignored.
     */
    static void checkStoichiometryNameMatching(
        const std::string& name, std::map<std::string, int>& stoich, 
        const std::vector<Element>& elements);

private:

    std::string m_name;
    
    std::string m_nasa7_name;
    std::string m_nasa9_name;
    
    ParticleRRHO* mp_rrho_model;
    
    double m_mw;
    int   m_charge;
    
    PhaseType m_phase;
    
    std::map<std::string, int> m_stoichiometry;
    
}; // class Species


/**
 * Swaps the contents of one Species object with another.
 */
void swap(Species& s1, Species& s2);

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_SPECIES_H

