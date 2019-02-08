/**
 * @file ParticleRRaHO.cpp
 *
 * @brief Implementation of the ParticleRRaHO class.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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


#include "ParticleRRaHO.h"
#include "Utilities.h"
#include "Constants.h"

#include <cassert>
#include <iostream>

using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ParticleRRaHO::ParticleRRaHO(const IO::XmlElement& xml_element)
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
	    double elec_energy;
            double diss_energy;

	    int e = 0; // electronic level index
            
            for ( ; level_iter != iter->end(); ++level_iter, e++) {
                if (level_iter->tag() == "level") {
                    level_iter->getAttribute("degeneracy", degeneracy);
                    level_iter->getAttribute("energy", temperature);

		    level_iter->getAttribute("elec_energy", elec_energy);
                    level_iter->getAttribute("diss_energy", diss_energy);

                    m_we.push_back(level_iter->getAttribute("we", we));
                    m_wexe.push_back(level_iter->getAttribute("wexe", wexe));
                    m_weye.push_back(level_iter->getAttribute("weye", weye));
                    m_weze.push_back(level_iter->getAttribute("weze", weze));
                    m_be.push_back(level_iter->getAttribute("Be", be));
                    m_ae.push_back(level_iter->getAttribute("Ae", ae));
                    m_de.push_back(level_iter->getAttribute("De", de));
                    
                    // convert from 1/cm to K
                    //m_electronic_energies.push_back(
                    //    std::make_pair(degeneracy, temperature * 1.4387));

		    m_electronic_energies.push_back(std::make_pair(degeneracy,
                        elec_energy));

		    m_electronical_energies.push_back(elec_energy);
		    m_dissociation_energies.push_back(diss_energy);
		    m_vibrational_temperatures.push_back(HP * C0 * m_we.at(e) / KB);

		    // Vibrational energy (using Morse potential)
		    double morse = 0.0;
                    int v = 0; // vibrational level index

		    vib_energy.push_back(std::vector<double>());
                    while (morse < (m_dissociation_energies.at(e) -
                                    m_electronical_energies.at(e)))
                    {
                        // Morse potential, Eq. 1.3, Nagnibeda & Kustova 2009.
                        morse = HP * C0 * (m_we.at(e) * (v + 0.5) -
                            m_wexe.at(e) * (v + 0.5) * (v + 0.5) +
                            m_weye.at(e) * (v + 0.5) * (v + 0.5) * (v + 0.5) +
                            m_weze.at(e) * (v + 0.5) * (v + 0.5) * (v + 0.5) * (v + 0.5));

                        m_vibrational_energies.push_back(morse);
                        vib_energy[e].push_back(morse);
                        v++;
                    }
		    m_vibrational_levels.push_back(v);

		    // Rotational energy
		    double rot = 0.0;
                    int r = 0; // rotational level index

		    while (rot < (m_dissociation_energies.at(e) -
                                  m_electronical_energies.at(e)))
                    {
                        rot = HP * C0 * m_be.at(e) * r * (r+1);
                        r++;
                    }
		    m_rotational_levels.push_back(r);
		    rot_energy.push_back(std::vector<std::vector<double>>());
                    for (v=0; v!=m_vibrational_levels.at(e); ++v) {
                        rot_energy[e].push_back(std::vector<double>());
                        for (r=0; r!=m_rotational_levels.at(e); ++r) {
                            rot = HP * C0 * (m_be.at(e) * r * (r+1));
                            rot_energy[e][v].push_back(rot);
                            m_rotational_energies.push_back(rot);
                        }
                    }
                }
            }
        }
	else if (iter->tag() == "mass") {
	    std::string str;
            iter->getAttribute("units", str, "kg");
            m_mass = Units(str).convertToBase(
            atof(iter->text().c_str()));
            //cout << " Particle's mass is: " << m_mass << endl;
	}
	else if (iter->tag() == "diameter") {
            m_diameter = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's diameter is: " << m_diameter << endl;
        }
        else if (iter->tag() == "LJeps") {
            m_eps_LJ = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's eps_LJ is: " << m_eps_LJ << endl;
        }
	else if (iter->tag() == "formation_energy") {
            m_formation_energy = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's formation energy is: " << m_formation_energy << endl;
        }
	else if (iter->tag() == "ionization_potential") {
             m_ionization_potential = atof(String::trim(iter->text()).c_str());
             //cout << " Particle's ionization potential is: "
             //   << m_ionization_potential << endl;
        }
        else if (iter->tag() == "inertia_moment") {
            m_inertia_moment = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's moment of inertia is: " << m_inertia_moment << endl;
        }
	else if (iter->tag() == "parker_Zeta_inf") {
            m_parker_zeta_inf = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's Parker zeta inf is: " << m_parker_zeta_inf << endl;
        }
	else if (iter->tag() == "reduced_oscillator_mass") {
            m_reduced_oscillator_mass = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's reduced oscillator mass is: " << m_reduced_oscillator_mass << endl;
        }
        else if (iter->tag() == "internuclear_distance") {
            m_r_e = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's internuclear distance is: " << m_r_e << endl;
        }
        else if (iter->tag() == "mA_mAB") {
            m_mA_mAB_ratio = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's mA/mAB is: " << m_mA_mAB_ratio << endl;
        }
	else if (iter->tag() == "mB_mAB") {
            m_mB_mAB_ratio = atof(String::trim(iter->text()).c_str());
            //cout << " Particle's mB/mAB is: " << m_mB_mAB_ratio << endl;
        }
    }
}

//==============================================================================

ParticleRRaHO::ParticleRRaHO(const ParticleRRaHO& rraho, const size_t level)
    : m_hform(rraho.m_hform + rraho.electronicEnergy(level).second),
      m_steric(rraho.m_steric),
      m_linearity(rraho.m_linearity),
      m_rotational_t(rraho.m_rotational_t),
      m_electronic_energies(1, rraho.electronicEnergy(level)),
      m_vibrational_energies(rraho.m_vibrational_energies)
{
    // Make sure the level used is actually present in the given RRaHO parameters
    assert(level < rraho.nElectronicLevels());
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation

