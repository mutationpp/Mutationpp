
#include "ParticleRRHO.h"
#include "Utilities.h"
#include "Constants.h"

#include <cassert>
#include <iostream>

using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ParticleRRHO::ParticleRRHO(const IO::XmlElement& xml_element)
    : m_hform(0), m_steric(0), m_linearity(0), m_rotational_t(0)
{
    // Load information stored in child elements
    IO::XmlElement::const_iterator iter = xml_element.begin();
    
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
            IO::XmlElement::const_iterator level_iter = iter->begin();
            
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

ParticleRRHO::ParticleRRHO(const ParticleRRHO& rrho, const size_t level)
    : m_hform(rrho.m_hform + rrho.electronicEnergy(level).second),
      m_steric(rrho.m_steric),
      m_linearity(rrho.m_linearity),
      m_rotational_t(rrho.m_rotational_t),
      m_electronic_energies(1, rrho.electronicEnergy(level)),
      m_vibrational_energies(rrho.m_vibrational_energies)
{
    // Make sure the level used is actually present in the given RRHO parameters
    assert(level < rrho.nElectronicLevels());
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation

