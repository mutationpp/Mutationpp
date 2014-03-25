
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
    std::swap(opt1.m_composition_setter, opt2.m_composition_setter);
    std::swap(opt1.m_has_default_composition, opt2.m_has_default_composition);
    std::swap(opt1.m_source, opt2.m_source);
    std::swap(opt1.m_state_model, opt2.m_state_model);
    std::swap(opt1.m_thermo_db, opt2.m_thermo_db);
    std::swap(opt1.m_mechanism, opt2.m_mechanism);
    std::swap(opt1.m_viscosity, opt2.m_viscosity);
    std::swap(opt1.m_thermal_conductivity, opt2.m_thermal_conductivity);
}

MixtureOptions::MixtureOptions()
    : m_composition_setter(m_default_composition),
      m_has_default_composition(false),
      m_source()
{ 
    setDefaultOptions();
}

MixtureOptions::MixtureOptions(const std::string& mixture)
    : m_composition_setter(m_default_composition),
      m_has_default_composition(false)
{
    loadFromFile(mixture);
}

MixtureOptions::MixtureOptions(const char* mixture)
    : m_composition_setter(m_default_composition),
      m_has_default_composition(false)
{
    loadFromFile(string(mixture));
}

MixtureOptions::MixtureOptions(IO::XmlElement& element)
    : m_composition_setter(m_default_composition),
      m_has_default_composition(false)
{
    loadFromXmlElement(element);
}


void MixtureOptions::setDefaultOptions()
{
    m_species_descriptor = "";
    m_source = "";
    m_state_model = "EquilTP";
    m_thermo_db   = "RRHO";
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
    ifstream file(m_source.c_str(), ios::in);

    // If that doesn't work then assume it is in MPP_DATA_DIRECTORY/mixtures.
    if (!file.is_open())
		m_source = getEnvironmentVariable("MPP_DATA_DIRECTORY") +
			"/mixtures/" + m_source;
    file.close();

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

