/**
 * @file SpeciesListDescriptor.cpp
 *
 * @brief Implementation of the SpeciesListDescriptor class.
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

#include "SpeciesListDescriptor.h"
#include "Species.h"
#include "Utilities.h"

#include <iostream>
#include <iterator>
using namespace std;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

SpeciesListDescriptor::SpeciesListDescriptor(std::string descriptor)
    : m_gases(false), m_solids(false), m_liquids(false)
{
    // First step is to look for implicit species defined by "{ rules }"
    size_t is = descriptor.find_first_of('{');
    if (is != std::string::npos) {
        size_t es = descriptor.find_first_of('}', is);
        std::string rule = descriptor.substr(is+1,es-is-1);
        descriptor.erase(is, es-is+1);

        // Split the rules up
        std::vector<std::string> rules;
        Utilities::String::tokenize(rule, rules, "with", false);
        
        // First check phase rule
        std::vector<std::string> phases;
        Utilities::String::tokenize(rules[0], phases, ", ");
        
        for (size_t i = 0; i < phases.size(); ++i) {
            if (phases[i] == "gases")
                m_gases = true;
            else if (phases[i] == "liquids")
                m_liquids = true;
            else if (phases[i] == "solids")
                m_solids = true;
            else if (phases[i] == "condensed") {
                m_solids = true;
                m_liquids = true;
            } else if (phases[i] == "all") {
                m_gases = true;
                m_solids = true;
                m_liquids = true;
            } else {
                throw InvalidInputError("species descriptor", descriptor)
                    << "Unknown phase keyword in implicit species rule.  "
                    << "Possible phase descriptors are 'gases', 'liquids', "
                    << "'solids', 'condensed', and 'all'.";
            }
        }
        
        // Second the elements
        std::vector<std::string> elements;
        Utilities::String::tokenize(rules[1], elements, ", ");
        m_element_names.insert(elements.begin(), elements.end());
    }
    
    // Separate out the species names
    separateSpeciesNames(descriptor);
    
    // Check on any species whose excited electronic states should be expanded
    for (int i = 0; i < m_species_names.size(); ++i) {
        const size_t size = m_species_names[i].size();
        if (size > 3 && m_species_names[i].substr(size-3,3) == "(*)") {
            m_species_names[i] = m_species_names[i].substr(0, size-3);
            m_expand_states.insert(m_species_names[i]);
        }
    }
}

//==============================================================================

void SpeciesListDescriptor::separateSpeciesNames(std::string descriptor)
{
    // First trim the descriptor string
    descriptor = Utilities::String::trim(descriptor);

    // State-machine with 2 states (either in quotes or not)
    bool in_quotes = (descriptor[0] == '\"');
    std::string name = "";

    for (int i = (in_quotes ? 1 : 0); i < descriptor.length(); ++i) {
        char c = descriptor[i];

        if (in_quotes) {
            // Add everything to the name until we escape out of the quotes
            switch(c) {
            case '\"':
                // Escape and add name to the list
                in_quotes = false;
                if (name.length() > 0) {
                    m_species_names.push_back(name);
                    name = "";
                }
                break;
            default:
                // Otherwise just add the character to the name
                name += c;
            }
        } else {
            switch(c) {
            // All white-space should be ignored
            case ' ': case '\n': case '\r': case '\t': case '\f': case '\v':
                if (name.length() > 0) {
                    m_species_names.push_back(name);
                    name = "";
                }
                break;
            // Cannot include quotation mark in a name
            case '\"':
                if (name.length() > 0) {
                    throw InvalidInputError("species name", name)
                        << "Cannot include quotation mark in species name.\n"
                        << "    " << descriptor.substr(0, i+1) << " <--";
                }
                in_quotes = true;
                break;
            default:
                name += c;
            }
        }
    }

    // Push back the last name
    if (name.length() > 0)
        m_species_names.push_back(name);
}

//==============================================================================

bool SpeciesListDescriptor::matches(const Species& species) const
{
    // Check if this species is present in the explicit list (includes excited
    // state species which were implicitly defined using the '*' character)
    for (int i = 0; i < m_species_names.size(); ++i) {
        const std::string& name = m_species_names[i];
        bool expand = (m_expand_states.count(name) == 0 ? false : true);
        
        if (species.name() == name)
            return !expand;
        
        if (species.groundStateName() == name)
            return expand;
    }
    
    // Check if the species is described by an implicit rule
    // Check phase
    if (species.phase() == GAS    && !m_gases)   return false;
    if (species.phase() == SOLID  && !m_solids)  return false;
    if (species.phase() == LIQUID && !m_liquids) return false;
    
    // Check elements
    Species::StoichList::const_iterator iter =
        species.stoichiometry().begin();
    for ( ; iter != species.stoichiometry().end(); ++iter)
        if (m_element_names.count(iter->first) == 0) return false;
    
    return true;
}

//==============================================================================

/// Predicate returns true if species name equals given name.
struct NameEquals {
    NameEquals(const std::string& str) : name(str) { }
    bool operator()(const Species& s) { return s.name() == name; }
    std::string name;
};

/// Predicate returns true if species ground state name equals given name.
struct GroundStateNameEquals {
    GroundStateNameEquals(const std::string& str) : name(str) { }
    bool operator()(const Species& s) { return s.groundStateName() == name; }
    std::string name;
};

/// Predicate returns true if species type equals given type.
struct TypeEquals {
    TypeEquals(ParticleType t) : type(t) { }
    bool operator()(const Species& s) { return s.type() == type; }
    ParticleType type;
};

/// Predicate returns true if species is condensed.
struct IsCondensed {
    bool operator()(const Species& s) { return s.phase() != GAS; }
};

void SpeciesListDescriptor::order(
    std::list<Species>& input, std::vector<Species>& output,
    std::vector<std::string>& missing) const
{
    std::list<Species> ordered_list;
    
    int index = 0;
    std::list<Species>::iterator it;
    
    // First order all of the species that are explicitly listed
    for (std::size_t i = 0; i < m_species_names.size(); ++i) {
        const std::string& name = m_species_names[i];
        
        // Don't expand name
        if (m_expand_states.find(name) == m_expand_states.end()) {
            it = std::find_if(input.begin(), input.end(), NameEquals(name));
            
            if (it == input.end())
                missing.push_back(name);
            else {
                ordered_list.push_back(*it);
                input.erase(it);
            }
        // Expand name
        } else {
            it = std::find_if(
                input.begin(), input.end(), GroundStateNameEquals(name));
            
            std::list<Species>::iterator lowest_level = it;
            while (lowest_level != input.end()) {
                it++;
                it = std::find_if(it, input.end(), GroundStateNameEquals(name));
                while (it != input.end()) {
                    if (it->level() < lowest_level->level())
                        lowest_level = it;
                    it++;
                    it = std::find_if(
                        it, input.end(), GroundStateNameEquals(name));
                }
            
                ordered_list.push_back(*lowest_level);
                input.erase(lowest_level);
                
                it = std::find_if(
                    input.begin(), input.end(), GroundStateNameEquals(name));
                lowest_level = it;
            }
        }
    }
    
    // Return early if we are missing some species as this is an error
    if (missing.size() > 0)
        return;
    
    // Now all that remains are the species that were implicitly defined
    ordered_list.insert(ordered_list.end(), input.begin(), input.end());
    input.clear();
    
    // Place the electron first
    it = find_if(
         ordered_list.begin(), ordered_list.end(), TypeEquals(ELECTRON));

    if (it != ordered_list.begin() && it != ordered_list.end())
        ordered_list.splice(ordered_list.begin(), ordered_list, it);

    // Move all condensed species to the end of the list
    int ncond = count_if(
        ordered_list.begin(), ordered_list.end(), IsCondensed());

    while (ncond-- > 0) {
        it = std::find_if(
            ordered_list.begin(), ordered_list.end(), IsCondensed());
        ordered_list.splice(ordered_list.end(), ordered_list, it);
    }

    // Finally copy the list to the output vector
    output.assign(ordered_list.begin(), ordered_list.end());
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
