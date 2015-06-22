/**
 * @file Species.h
 *
 * @brief Declaration of the Species class.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef THERMO_SPECIES_H
#define THERMO_SPECIES_H

#include <map>
#include <string>
//#include <set>
#include <vector>

// Forward declaration of the XmlElement Class
namespace Mutation {
namespace Utilities {
namespace IO {
    class XmlElement;
} } }
     

namespace Mutation {
    namespace Thermodynamics {
    
//class ParticleRRHO;

/**
 * Responsible for loading element data from the elements.xml file.
 */
class Element
{
public:
    
    /**
     * Loads an element from an "element" XmlElement node.
     */
    Element(const Mutation::Utilities::IO::XmlElement& xml_element);
    
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
    double      m_atomic_mass;
    int         m_charge;

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
 * Stores basic information about a species including how many elements belong
 * to the species, its charge, and phase type.
 */
class Species
{
public:

    typedef std::vector< std::pair<std::string, int> > StoichList;
    
    Species()
        : m_name(""), m_ground_state_name(""), m_mw(0.0), m_charge(0),
          m_phase(GAS), m_type(ATOM), m_level(0)
    { }
    
    
    /**
     * Create a new species object.
     */
    Species(
        const std::string& name,
        const PhaseType phase,
        const StoichList& stoichiometry,
        const std::vector<Element>& elements);
    
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
    ~Species(){};
    
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
     * Returns the name of the ground state for this species.
     */
    const std::string& groundStateName() const {
        return m_ground_state_name;
    }
    
    /**
     * Returns type of species (ie: electron, atom, or molecule).
     *
     * @see enum ParticleType
     */
    ParticleType type() const {
        return m_type;
    }
    
    /**
     * Returns the molecular weight of the species in kg/mol.
     */
    double molecularWeight() const {
        return m_mw;
    }
    
    /**
     * Returns the number of atoms of a particular element contained in this
     * species.  For example calling nAtoms("C") for the C2H2 molecule would
     * return 2.
     */
    int nAtoms(const std::string &element) const {
        StoichList::const_iterator iter = m_stoichiometry.begin();
        for ( ; iter != m_stoichiometry.end(); ++iter)
            if (iter->first == element) return iter->second;
        return 0;
    }
    
    /**
     * Returns the total number of atoms in this species.  For example, C2H2 has
     * 4 total atoms (2 Carbon and 2 Hydrogen atoms).  Note that this does not
     * include electrons in the calculation.
     */
    int nAtoms() const {
        int atoms = 0;
        StoichList::const_iterator iter = m_stoichiometry.begin();
        for ( ; iter != m_stoichiometry.end(); ++iter)
            atoms += (iter->first != "e-" ? iter->second : 0);
        return atoms;
    }
    
    /**
     * Returns a vector of element name/number pairs for the elements that make
     * up this species.
     */
    const StoichList& stoichiometry() const {
        return m_stoichiometry;
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
     * Returns the excited state level if this is an excited state species.
     */
    std::size_t level() const {
        return m_level;
    }
    
    friend void swap(Species&, Species&);

private:

    /**
     * Determines the stoichiometry of this species from a stoichiometry string.
     */
    //void loadStoichiometry(
    //    const std::string& stoichiometry, const std::vector<Element> &elements);

    /**
     * Will check the species name against its stoichiometry to determine if
     * the stoichiometry given by the user matches the stoichiometry found by
     * parsing the name.  Note that the name can be anything, therefore if it is
     * not in the standard form then it will simply be ignored.
     */
    static void checkStoichiometryNameMatching(
        const std::string& name, std::map<std::string, int>& stoich, 
        const std::vector<Element>& elements);

private:

    std::string  m_name;
    std::string  m_ground_state_name;
    double       m_mw;
    int          m_charge;
    PhaseType    m_phase;
    ParticleType m_type;
    std::size_t  m_level;
    
    StoichList m_stoichiometry;
    
    //ParticleRRHO* mp_rrho_model;
    
}; // class Species


/**
 * Swaps the contents of one Species object with another.
 */
void swap(Species& s1, Species& s2);

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_SPECIES_H

