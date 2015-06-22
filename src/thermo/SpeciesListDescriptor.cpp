/**
 * @file SpeciesListDescriptor.cpp
 *
 * @brief Implementation of the SpeciesListDescriptor class.
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

#include "SpeciesListDescriptor.h"
#include "Species.h"
#include "Utilities.h"

#include <iostream>
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
                std::cout << "Unknown phase keyword in implicit species rule"
                          << std::endl;
                exit(1);
            }
        }
        
        // Second the elements
        std::vector<std::string> elements;
        Utilities::String::tokenize(rules[1], elements, ", ");
        m_element_names.insert(elements.begin(), elements.end());
    }
    
    // Tokenize the species list to get the species names
    Utilities::String::tokenize(
        Utilities::String::trim(descriptor), m_species_names, " \n\r\t\f\v");
    
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

void SpeciesListDescriptor::order(
    std::list<Species>& input, std::vector<Species>& output,
    std::vector<std::string>& missing) const
{
    // Prepare memory for the ordered list
    output.resize(input.size());
    
    int index = 0;
    std::list<Species>::iterator species_iter;
    
    // First order all of the species that are explicitly listed
    for (std::size_t i = 0; i < m_species_names.size(); ++i) {
        const std::string& name = m_species_names[i];
        
        if (m_expand_states.count(name) == 0) {
            species_iter = input.begin();
            while (species_iter != input.end() && species_iter->name() != name)
                species_iter++;
            
            if (species_iter == input.end())
                missing.push_back(name);
            else {
                output[index++] = *species_iter;
                input.erase(species_iter);
            }
        } else {
            species_iter = input.begin();
            while (species_iter->groundStateName() != name &&
                   species_iter != input.end())
                species_iter++;
            
            std::list<Species>::iterator lowest_level = species_iter;
            while (lowest_level != input.end()) {
                while (species_iter != input.end()) {
                    if (species_iter->groundStateName() == name &&
                        species_iter->level() < lowest_level->level())
                        lowest_level = species_iter;
                    species_iter++;
                }
            
                output[index++] = *lowest_level;
                input.erase(lowest_level);
                
                species_iter = input.begin();
                while (species_iter->groundStateName() != name &&
                       species_iter != input.end())
                    species_iter++;
                lowest_level = species_iter;
            }
        }
    }
    
    // Return early if we are missing some species as this is an error
    if (missing.size() > 0)
        return;
    
    // Now all that remains are the species that were implicitly defined
    species_iter = input.begin();
    while (species_iter != input.end()) {
        output[index++] = *species_iter;
        species_iter = input.erase(species_iter);
    }
    
    // Place the electron first
    for (index = 0; index < output.size(); ++index)
        if (output[index].type() == ELECTRON) break;

    if (index < output.size() && index != 0) {
        const Species electron = output[index];
        while (index > 0) {
            output[index] = output[index-1];
            index--;
        }
        output[0] = electron;
    }

    // Move all condensed species to the end of the list
    int end;
    index = end = output.size()-1;

    while (index > -1) {
        if (output[index].phase() != GAS)
            swap(output[index], output[end--]);
        index--;
    }
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
