
#include "ParticleRRHO.h"
#include "Utilities.h"

#include <cassert>
#include <iostream>

using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ParticleRRHO::ParticleRRHO(IO::XmlElement& xml_element)
    : m_hform(0), m_steric(0), m_linearity(0), m_rotational_t(0)
{
    // Load information stored in child elements
    IO::XmlElement::Iterator iter = xml_element.begin();        
    
    for ( ; iter != xml_element.end(); ++iter) {
        if (iter->tag() == "formation_enthalpy") {
            m_hform = atof(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "steric_factor") {
            m_steric = atoi(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "linear") {
            std::string yesno = String::trim(iter->text());
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
            m_rotational_t = atof(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "vibrational_temperatures") {
            std::vector<std::string> tokens;
            String::tokenize(iter->text(), tokens, " ,\t\n\r\f\v");
            
            std::vector<std::string>::const_iterator t_iter = tokens.begin();
            for ( ; t_iter != tokens.end(); ++t_iter)
                m_vibrational_energies.push_back(atof(t_iter->c_str()));
        }
        else if (iter->tag() == "electronic_levels") {
            IO::XmlElement::Iterator level_iter = iter->begin();
            
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

//==============================================================================

ParticleRRHO::ParticleRRHO(const ParticleRRHO* p_rrho, const size_t level)
    : m_hform(p_rrho->m_hform),
      m_steric(p_rrho->m_steric),
      m_linearity(p_rrho->m_linearity),
      m_rotational_t(p_rrho->m_rotational_t),
      m_electronic_energies(),
      m_vibrational_energies(p_rrho->m_vibrational_energies)
{
    // Make sure the level used is actually present in the given RRHO parameters
    assert(level < p_rrho->nElectronicLevels());
    
    // Copy only the level given
    m_electronic_energies.push_back(p_rrho->electronicEnergy(level));
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation

