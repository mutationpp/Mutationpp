/**
 * @file MixtureOptions.cpp
 *
 * @brief Provides MixtureOptions class definition.
 * @see Mutation::MixtureOptions
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

#include <iostream>
#include <fstream>

#include "MixtureOptions.h"
#include "Utilities.h"

using namespace std;
using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

namespace Mutation {

void swap(MixtureOptions& opt1, MixtureOptions& opt2)
{
    std::swap(opt1.m_species_descriptor, opt2.m_species_descriptor);
    std::swap(opt1.m_compositions, opt2.m_compositions);
    std::swap(opt1.m_default_composition, opt2.m_default_composition);
    //std::swap(opt1.m_has_default_composition, opt2.m_has_default_composition);
    std::swap(opt1.m_source, opt2.m_source);
    std::swap(opt1.m_state_model, opt2.m_state_model);
    std::swap(opt1.m_thermo_db, opt2.m_thermo_db);
    std::swap(opt1.m_mechanism, opt2.m_mechanism);
    std::swap(opt1.m_viscosity, opt2.m_viscosity);
    std::swap(opt1.m_thermal_conductivity, opt2.m_thermal_conductivity);
    std::swap(opt1.m_gsi_mechanism, opt2.m_gsi_mechanism);
}

MixtureOptions::MixtureOptions()
    : m_source()
{ 
    setDefaultOptions();
}

MixtureOptions::MixtureOptions(const std::string& mixture)
{
    loadFromFile(mixture);
}

MixtureOptions::MixtureOptions(const char* mixture)
{
    loadFromFile(string(mixture));
}

MixtureOptions::MixtureOptions(IO::XmlElement& element)
{
    loadFromXmlElement(element);
}


void MixtureOptions::setDefaultOptions()
{
    m_species_descriptor = "";
    m_default_composition = -1;
    m_source = "";
    m_state_model = "ChemNonEq1T";
    m_thermo_db   = "RRHO";
    m_mechanism   = "none";
    m_viscosity   = "Chapmann-Enskog_LDLT";
    m_thermal_conductivity = "Chapmann-Enskog_LDLT";
    m_gsi_mechanism = "none";
}

void MixtureOptions::loadFromFile(const string& mixture)
{
    // Initalize to the default options
    setDefaultOptions();
    
    // Get the mixture path on this computer.
    m_source = databaseFileName(mixture, "mixtures");

    // Now load the XML file
    IO::XmlDocument mixture_doc(m_source);
    IO::XmlElement root = mixture_doc.root();
    
    // Now we can load from the mixture XmlElement
    loadFromXmlElement(root);
}

void MixtureOptions::loadFromXmlElement(IO::XmlElement& element)
{
    if (element.tag() != "mixture")
        element.parseError(
            "XmlElement is not of 'mixture' type!");

    // Get the name of the mixture reaction mechanism
    element.getAttribute("mechanism", m_mechanism, m_mechanism);
    
    // Get the type of thermodynamic database to use
    element.getAttribute("thermo_db", m_thermo_db, m_thermo_db);
    
    // Get the viscosity algorithm
    element.getAttribute("viscosity", m_viscosity, m_viscosity);
    
    // Get the thermal conductivity algorithm
    element.getAttribute(
        "thermal_conductivity", m_thermal_conductivity, m_thermal_conductivity);
    
    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute("gsi_mechanism", m_gsi_mechanism, m_gsi_mechanism);

    // Get the state model
    element.getAttribute("state_model", m_state_model, m_state_model);

    // Loop over all of the mixture child elements
    IO::XmlElement::const_iterator iter;
    for (iter = element.begin(); iter != element.end(); ++iter) {
        // Load the species list
        if (iter->tag() == "species")
            m_species_descriptor = String::trim(iter->text());
        else if (iter->tag() == "element_compositions")
            loadElementCompositions(*iter);
    }
}

void MixtureOptions::loadElementCompositions(const IO::XmlElement& element)
{
    // Load the compositions
    IO::XmlElement::const_iterator comps = element.begin();
    while (comps != element.end()) {
        if (comps->tag() == "composition")
            if (!addComposition(Composition(*comps)))
                comps->parseError("Redefinition of a named composition!");
        comps++;
    }

    // Set the default composition
    std::string name = "";
    element.getAttribute("default", name, name);

    if (name != "") {
        if (!setDefaultComposition(name))
            element.parseError(
                "Default composition does not match one of named compositions!");
    }
}

bool MixtureOptions::addComposition(const Composition& c, bool make_default)
{
    // Check that this composition has a unique name
    const int n = m_compositions.size();
    for (int i = 0; i < n; ++i)
        if (c.name() == m_compositions[i].name())
            return false;

    m_compositions.push_back(c);
    if (make_default)
        m_default_composition = n;

    return true;
}

bool MixtureOptions::setDefaultComposition(const std::string& name)
{
    m_default_composition = -1;
    for (int i = 0; i < m_compositions.size(); ++i) {
        if (name == m_compositions[i].name()) {
            m_default_composition = i;
            break;
        }
    }

    return hasDefaultComposition();
}


} // namespace Mutation

