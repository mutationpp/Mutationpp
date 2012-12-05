
#include "ParticleRRHO.h"
#include "utilities.h"
#include "XMLite.h"

#include <iostream>

using namespace utils;

ParticleRRHO::ParticleRRHO(XmlElement& xml_element)
    : m_hform(0), m_steric(0), m_linearity(0), m_rotational_t(0)
{
    // Load information stored in child elements
    XmlElement::Iterator iter = xml_element.begin();        
    
    for ( ; iter != xml_element.end(); ++iter) {
        if (iter->tag() == "formation_enthalpy") {
            m_hform = atof(StringUtils::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "steric_factor") {
            m_steric = atoi(StringUtils::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "linear") {
            std::string yesno = StringUtils::trim(iter->text());
            if (yesno == "yes")
                m_linearity = 2;
            else if (yesno == "no")
                m_linearity = 3;
            else {
                iter->parseError(
                    "Values for linear can only be \"yes\" or \"no\"!");
            }
        }
        else if (iter->tag() == "rotational_temperature") {
            m_rotational_t = atof(StringUtils::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "vibrational_temperatures") {
            std::vector<std::string> tokens;
            StringUtils::tokenize(iter->text(), tokens, " ,\t\n\r\f\v");
            
            std::vector<std::string>::const_iterator t_iter = tokens.begin();
            for ( ; t_iter != tokens.end(); ++t_iter)
                m_vibrational_energies.push_back(atof(t_iter->c_str()));
        }
        else if (iter->tag() == "electronic_levels") {
            XmlElement::Iterator level_iter = iter->begin();
            
            int    degeneracy;
            double temperature;
            
            for ( ; level_iter != iter->end(); ++level_iter) {
                if (level_iter->tag() == "level") {
                    level_iter->getAttribute("degeneracy",  degeneracy);
                    level_iter->getAttribute("energy", temperature);
                    
                    // convert from 1/cm to K
                    m_electronic_energies.push_back(
                        std::make_pair(degeneracy, temperature * 1.4387));
                }
            }
        }
    }
}
