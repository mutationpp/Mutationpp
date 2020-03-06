/**
 * @file MillikanWhite.cpp
 *
 * @brief Implementation of classes related to Millikan and White model.
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


#include "MillikanWhite.h"
#include "Utilities.h"
#include "ParticleRRHO.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Mutation::Thermodynamics;

namespace Mutation {
    namespace Transfer {

//==============================================================================

MillikanWhiteVibrator::MillikanWhiteVibrator(
    const XmlElement& node, const class Thermodynamics& thermo)
{
    assert(node.tag() == "vibrator");
        
    // Get the name of this species
    std::string name;
    node.getAttribute("species", name, "must provide species name!");
    m_index  = thermo.speciesIndex(name);
    m_thetav = loadThetaV(name);
    
    // Get the limiting cross-section if available
    node.getAttribute("omegav", m_omegav, 3.0E-21);
    
    const Species& vibrator = thermo.species(name);
    
    // Loop over each heavy species in thermo
    int offset = (thermo.hasElectrons() ? 1 : 0);
    XmlElement::const_iterator partner_iter;
    
    double a, b, mu;
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        // Get collision partner
        const Species& partner = thermo.species(i+offset);
        
        // Compute reduced mass of this pair
        mu = (vibrator.molecularWeight() * partner.molecularWeight()) /
             (vibrator.molecularWeight() + partner.molecularWeight());
            
        // Use a and b data from data file or use the defaults if the pair
        // is not present in the file
        if ((partner_iter = node.findTagWithAttribute(
            "partner", "species", partner.name())) != node.end()) {
            // Get a, b from parnter node
            partner_iter->getAttribute("a", a, "must provide constant a!");
            partner_iter->getAttribute("b", b, "must provide constant b!");
            
            // Add Millikan-White data for collision pair
            m_partners.push_back(MillikanWhitePartner(a, b, mu));
        } else {
            if (m_thetav < 0.0) {
                std::cout << "Error: Did not find vibrational temperature for "
                          << name << "." << std::endl;
                std::exit(1);
            }

            // Add Millikan-White data using defaults
            m_partners.push_back(MillikanWhitePartner(mu, m_thetav));
        }
    }
}

//==============================================================================

MillikanWhiteVibrator::MillikanWhiteVibrator(
    const std::string& name, const class Thermodynamics& thermo)
    : m_omegav(3.0E-21),
      m_index(thermo.speciesIndex(name)),
      m_thetav(loadThetaV(name))
{
    // Rely on thetav for all partners so make sure it was found
    if (m_thetav < 0.0) {
        std::cout << "Error: Did not find vibrational temperature for "
                  << name << "." << std::endl;
        std::exit(1);
    }

    const Species& vibrator = thermo.species(name);
    
    // Loop over each heavy species in thermo
    int offset = (thermo.hasElectrons() ? 1 : 0);
    double mu;
    
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        // Get collision partner
        const Species& partner = thermo.species(i+offset);
        
        // Compute reduced mass of this pair
        mu = (vibrator.molecularWeight() * partner.molecularWeight()) /
             (vibrator.molecularWeight() + partner.molecularWeight());
        
        // Add Millikan-White data using defaults
        m_partners.push_back(MillikanWhitePartner(mu, m_thetav));
    }
}

//==============================================================================

double MillikanWhiteVibrator::loadThetaV(const std::string& name)
{
    // Get the species.xml path on this computer.
    std::string filename = databaseFileName("species.xml", "thermo");

    // Now look for the species
    XmlDocument doc(filename);
    XmlElement::const_iterator iter = doc.root().findTagWithAttribute(
        "species", "name", name);
    if (iter == doc.root().end()) return -1.0;

    // Look for the thermodynamic data
    XmlElement::const_iterator thermo_iter = iter->findTagWithAttribute(
        "thermodynamics", "type", "RRHO");
    if (thermo_iter == iter->end()) return -1.0;

    // Look for the vibrational temperature data
    iter = thermo_iter->findTag("vibrational_temperatures");
    if (iter == thermo_iter->end()) return -1.0;

    std::stringstream ss(iter->text());
    double thetav; ss >> thetav;

    return thetav;
}

//==============================================================================

MillikanWhite::MillikanWhite(const class Thermodynamics& thermo)
{
    // Get the VT.xml file location.
    std::string filename = databaseFileName("VT.xml", "transfer");
    
    // Open the VT file as an XML document
    XmlDocument doc(filename);
    
    // Find the Millikan-White data
    XmlElement::const_iterator iter = doc.root().findTag("Millikan-White");
    if (iter == doc.root().end())
        doc.root().parseError("Could not find Millikan-White element.");
    const XmlElement& root = *iter;

    // Loop over all molecules and load the Millikan-White data associated with
    // them
    int offset = (thermo.hasElectrons() ? 1 : 0);
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        const Species& species = thermo.species(i+offset);
        
        // If this molecule can vibrate, add it to the list
        if (species.type() == MOLECULE) {
            // If vibrator is not in in VT.xml take the characteristic vibrational
            // temperature from species.xml
            if ((iter = root.findTagWithAttribute(
                "vibrator", "species", species.name())) != root.end())
                m_vibrators.push_back(
                    MillikanWhiteVibrator(*iter, thermo));
            else
                m_vibrators.push_back(
                    MillikanWhiteVibrator(species.name(), thermo));
        }
    }
}

//==============================================================================

    } // namespace Transfer
} // namespace Mutation
