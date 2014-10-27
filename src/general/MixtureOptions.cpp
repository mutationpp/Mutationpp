/**
 * @file MixtureOptions.cpp
 *
 * @brief Provides MixtureOptions class definition.
 * @see Mutation::MixtureOptions
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

#include <iostream>
#include <fstream>

#include "MixtureOptions.h"
#include "Utilities.h"

using namespace std;
using namespace Mutation::Utilities;

namespace Mutation {

void swap(MixtureOptions& opt1, MixtureOptions& opt2)
{
    std::swap(opt1.m_species_descriptor, opt2.m_species_descriptor);
    std::swap(opt1.m_has_default_composition, opt2.m_has_default_composition);
    std::swap(opt1.m_load_transport, opt2.m_load_transport);
    std::swap(opt1.m_source, opt2.m_source);
    std::swap(opt1.m_state_model, opt2.m_state_model);
    std::swap(opt1.m_thermo_db, opt2.m_thermo_db);
    std::swap(opt1.m_mechanism, opt2.m_mechanism);
    std::swap(opt1.m_viscosity, opt2.m_viscosity);
    std::swap(opt1.m_thermal_conductivity, opt2.m_thermal_conductivity);
}

MixtureOptions::MixtureOptions()
    : m_has_default_composition(false),
      m_source()
{ 
    setDefaultOptions();
}

MixtureOptions::MixtureOptions(const std::string& mixture)
    : m_has_default_composition(false)
{
    loadFromFile(mixture);
}

MixtureOptions::MixtureOptions(const char* mixture)
    : m_has_default_composition(false)
{
    loadFromFile(string(mixture));
}

MixtureOptions::MixtureOptions(IO::XmlElement& element)
    : m_has_default_composition(false)
{
    loadFromXmlElement(element);
}


void MixtureOptions::setDefaultOptions()
{
    m_species_descriptor = "";
    m_source = "";
    m_state_model = "ChemNonEq1T";
    m_thermo_db   = "RRHO";
    m_load_transport = true;
    m_mechanism   = "none";
    m_viscosity   = "LDLT";
    m_thermal_conductivity = "LDLT";
}

void MixtureOptions::loadFromFile(const string& mixture)
{
    // Initalize to the default options
    setDefaultOptions();
    
    // Get the mixture path on this computer (first assume the local directory)
    m_source = mixture + ".xml";
    {
        ifstream file(m_source.c_str(), ios::in);

        // If that doesn't work then assume it is in MPP_DATA_DIRECTORY/mixtures.
        if (!file.is_open())
            m_source = getEnvironmentVariable("MPP_DATA_DIRECTORY") +
                "/mixtures/" + m_source;
    }

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
    
    // Get the state model
    element.getAttribute("state_model", m_state_model, m_state_model);
    
    // Check if we should load the transport data at all
    element.getAttribute("use_transport", m_load_transport, m_load_transport);

    // Loop over all of the mixture child elements
    IO::XmlElement::const_iterator iter;
    for (iter = element.begin(); iter != element.end(); ++iter) {
        // Load the species list
        if (iter->tag() == "species") {
            m_species_descriptor = String::trim(iter->text());
        
        // Load the default element fractions, note that we only check for valid
        // format, not for valid elements or fractions, this is left up to the
        // class that uses this information
        } else if (iter->tag() == "default_element_fractions") {
            vector<string> element_strings;
            String::tokenize(iter->text(), element_strings, ":, \n\r\t\f");
            
            if (element_strings.size() % 2 != 0) {
                iter->parseError(
                    "Default element fractions should have the format:\n" +
                    string("   name : fraction, name : fraction, ..."));
            }
            
            m_default_composition.clear();
            for (int i = 0; i < element_strings.size(); i+=2) {
                if (String::isNumeric(element_strings[i+1])) {            
                    m_default_composition.push_back(
                        make_pair(element_strings[i], 
                        atof(element_strings[i+1].c_str())));
                } else {
                    iter->parseError(
                        "Element fraction should be a real value!");
                }
            }
            
            m_has_default_composition = true;
        }
    }
}


} // namespace Mutation

