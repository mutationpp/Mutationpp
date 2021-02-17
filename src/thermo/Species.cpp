/**
 * @file Species.cpp
 *
 * @brief Implementation of the Species class.
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

#include "Errors.h"
#include "Species.h"
#include "SpeciesNameFSM.h"
#include "Utilities.h"
#include <algorithm>

using namespace std;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

const std::vector<Element>& Element::database()
{
    static std::vector<Element> elements;

    // If we have already loaded then simply return the vector
    if (elements.size() > 0)
        return elements;

    // Load from the XML database
    IO::XmlDocument element_doc(databaseFileName("elements.xml", "thermo"));
    IO::XmlElement::const_iterator element_iter = element_doc.root().begin();

    for ( ; element_iter != element_doc.root().end(); ++element_iter)
        elements.push_back(*element_iter);

    return elements;
}

//==============================================================================

Element::Element(const IO::XmlElement &xml_element)
{
    xml_element.getAttribute("name", m_name, 
        "Element must have a name attribute!");
    
    xml_element.getAttribute("charge", m_charge, 0);
    
    IO::XmlElement::const_iterator child_iter = xml_element.begin();
    IO::XmlElement::const_iterator child_end  = xml_element.end();
    
    for ( ; child_iter != child_end; ++child_iter) {
        if (child_iter->tag() == "mw") {
            std::string str;
            child_iter->getAttribute("units", str, "kg/mol");
            m_atomic_mass = Units(str).convertToBase(
                atof(child_iter->text().c_str()));
        }
    }
}

//==============================================================================

Species::Species(const Species& to_copy)
    : m_name(to_copy.m_name),
      m_ground_state_name(to_copy.m_ground_state_name),
      m_mw(to_copy.m_mw),
      m_charge(to_copy.m_charge),
      m_phase(to_copy.m_phase),
      m_type(to_copy.m_type),
      m_level(to_copy.m_level),
      m_stoichiometry(to_copy.m_stoichiometry)
{ }

//==============================================================================

Species::Species(const Species& to_copy, const size_t level) :
    m_name(to_copy.m_name),
    m_ground_state_name(to_copy.m_name),
    m_mw(to_copy.m_mw),
    m_charge(to_copy.m_charge),
    m_phase(to_copy.m_phase),
    m_type(to_copy.m_type),
    m_level(level),
    m_stoichiometry(to_copy.m_stoichiometry)
{ 
    stringstream ss;
    ss << "(" << m_level << ")";
    m_name += ss.str();
}

//==============================================================================

Species::Species(
    const std::string& name, const PhaseType phase) :
    m_name(name),
    m_ground_state_name(name),
    m_mw(0.0),
    m_charge(0),
    m_phase(phase),
    m_type(ATOM),
    m_level(0)
{
    // Parse the stoichiometry from the name
    SpeciesNameFSM sm;
    if (!sm.parse(name)) {
        throw InvalidInputError ("species name", name)
            << "If the stoichiometry is not explicitly given, then species "
            << "name must be formatted as a valid molecular formula.";
    }

    m_stoichiometry.assign(
        sm.stoichiometry().begin(), sm.stoichiometry().end());

    // Initialize everything else from the stoichiometry
    initDataFromStoichiometry();
}

//==============================================================================

Species::Species(
    const std::string& name, const PhaseType phase,
    const StoichList& stoichiometry) :
    m_name(name),
    m_ground_state_name(name),
    m_mw(0.0),
    m_charge(0),
    m_phase(phase),
    m_type(ATOM),
    m_level(0),
    m_stoichiometry(stoichiometry)
{
    initDataFromStoichiometry();
}

//==============================================================================

Species::Species(const Mutation::Utilities::IO::XmlElement& xml_element) :
    m_mw(0.0),
    m_charge(0),
    m_type(ATOM),
    m_level(0)
{
    // Species name
    xml_element.getAttribute("name", m_name,
        "Species must have a name attribute!");
    m_ground_state_name = m_name;

    // Phase
    std::string phase_str;
    xml_element.getAttribute("phase", phase_str, string("gas"));
    phase_str = String::toLowerCase(phase_str);

    if (phase_str == "gas")
        m_phase = GAS;
    else if (phase_str == "liquid")
        m_phase = LIQUID;
    else if (phase_str == "solid")
        m_phase = SOLID;
    else
        xml_element.parseError(
            "Invalid phase description for species \"" + m_name +
            "\", must be \"gas\" (default), \"liquid\", or \"solid\".");

    // Elemental stoichiometry
    IO::XmlElement::const_iterator stoich_iter =
        xml_element.findTag("stoichiometry");
    std::string stoich_str = stoich_iter->text();

    // Split up the string in the format "A:a, B:b, ..." into a vector of
    // strings matching ["A", "a", "B", "b", ... ]
    vector<string> stoich_tokens;
    String::tokenize(stoich_str, stoich_tokens, " \t\f\v\n\r:,");

    // Check that the vector has an even number of tokens (otherwise there must
    // be an error in the syntax)
    if (stoich_tokens.size() % 2 != 0) {
        xml_element.parseError(
            "Error in species \"" + m_name +
            "\" stoichiometry definition, invalid syntax!");
    }

    for (int i = 0; i < stoich_tokens.size(); i+=2)
        m_stoichiometry.push_back(std::make_pair(
            stoich_tokens[i], atoi(stoich_tokens[i+1].c_str())));

    initDataFromStoichiometry();
}

//==============================================================================

void Species::initDataFromStoichiometry()
{
    m_mw = 0.0;
    m_charge = 0;

    StoichList::const_iterator iter = m_stoichiometry.begin();
    for ( ; iter != m_stoichiometry.end(); ++iter) {
        std::vector<Element>::const_iterator el = Element::database().begin();
        for ( ; el != Element::database().end(); ++el)
            if (el->name() == iter->first) break;

        if (el == Element::database().end()) {
            throw InvalidInputError("element name", iter->first)
                ("species name", m_name)
                << "This element is not provided in the element database, "
                << "therefore it cannot be used to define the species "
                << "stoichiometry.";
        }

        m_mw += iter->second * el->atomicMass();
        m_charge += iter->second * el->charge();
    }

    // Determine the species type
    int atoms = nAtoms();
    m_type = (atoms == 0 ? ELECTRON : (atoms == 1 ? ATOM : MOLECULE));
}

//==============================================================================

void swap(Species& s1, Species& s2)
{
    std::swap(s1.m_name, s2.m_name);
    std::swap(s1.m_ground_state_name, s2.m_ground_state_name);
    std::swap(s1.m_mw, s2.m_mw);
    std::swap(s1.m_charge, s2.m_charge);
    std::swap(s1.m_phase, s2.m_phase);
    std::swap(s1.m_type, s2.m_type);
    std::swap(s1.m_level, s2.m_level);
    std::swap(s1.m_stoichiometry, s2.m_stoichiometry);
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation


