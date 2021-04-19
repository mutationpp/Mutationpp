/**
 * @file FHO.cpp
 *
 * @brief Implementation of classes related to FHO model.
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

#include "FHO.h"
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

FHOVibrator::FHOVibrator(
    const XmlElement& node, const class Thermodynamics& thermo)
{
    assert(node.tag() == "vibrator");
        
    // Get the name of this species
    std::string name;
    node.getAttribute("species", name, "must provide species name!");
    //m_index  = thermo.speciesIndex(name);
    //m_thetav = loadThetaV(name);
    
    // Get the limiting cross-section if available
    //node.getAttribute("omegav", m_omegav, 3.0E-21);
    
    const Species& vibrator = thermo.species(name);
    
    // Loop over each heavy species in thermo
    int offset = (thermo.hasElectrons() ? 1 : 0);
    XmlElement::const_iterator partner_iter;
    
    double beta, E_Morse, svt;
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        // Get collision partner
        const Species& partner = thermo.species(i+offset);
        
        // Compute reduced mass of this pair
        //mu = (vibrator.molecularWeight() * partner.molecularWeight()) /
        //     (vibrator.molecularWeight() + partner.molecularWeight());
            
        // Use beta, E_Morse and svt data from data file or use the defaults 
	// if the pair is not present in the file
        if ((partner_iter = node.findTagWithAttribute(
            "partner", "species", partner.name())) != node.end()) {
            // Get abeta, E_Morse and svt from parnter node
            partner_iter->getAttribute("beta", beta, "must provide constant beta!");
            partner_iter->getAttribute("E_Morse", E_Morse, "must provide constant E_Morse!");
            partner_iter->getAttribute("svt", svt, "must provide constant svt!");
            
            // Add FHO data for collision pair
            m_partners.push_back(FHOPartner(beta, E_Morse, svt));
        } else {
            //if (m_thetav < 0.0) {
                std::cout << "Error: Did not find vibrational temperature for "
                          << name << "." << std::endl;
                std::exit(1);
            //}

            // Add FHO data using defaults
            m_partners.push_back(FHOPartner());
        }
    }
}

//==============================================================================

FHOVibrator::FHOVibrator(
    const std::string& name, const class Thermodynamics& thermo)
    //: m_omegav(3.0E-21),
    //  m_index(thermo.speciesIndex(name)),
    //  m_thetav(loadThetaV(name))
{
    // Rely on thetav for all partners so make sure it was found
    //if (m_thetav < 0.0) {
    //    std::cout << "Error: Did not find vibrational temperature for "
    //              << name << "." << std::endl;
    //    std::exit(1);
    //}

    const Species& vibrator = thermo.species(name);
    
    // Loop over each heavy species in thermo
    int offset = (thermo.hasElectrons() ? 1 : 0);
    //double mu;
    
    for (int i = 0; i < thermo.nHeavy(); ++i) {
        // Get collision partner
        const Species& partner = thermo.species(i+offset);
        
        // Compute reduced mass of this pair
        //mu = (vibrator.molecularWeight() * partner.molecularWeight()) /
        //     (vibrator.molecularWeight() + partner.molecularWeight());
        
        // Add FHO data using defaults
        m_partners.push_back(FHOPartner());
    }
}

//==============================================================================

FHO::FHO(const class Thermodynamics& thermo)
{
    // Get the VT_FHO.xml file location.
    std::string filename = databaseFileName("VT_FHO.xml", "transfer");
    
    // Open the VT_FHO file as an XML document
    XmlDocument doc(filename);
    
    // Find the FHO data
    XmlElement::const_iterator iter = doc.root().findTag("FHO");
    if (iter == doc.root().end())
        doc.root().parseError("Could not find FHO element.");
    const XmlElement& root = *iter;

    // Loop over all molecules and load the FHO data associated with them
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
                    FHOVibrator(*iter, thermo));
            else
                m_vibrators.push_back(
                    FHOVibrator(species.name(), thermo));
        }
    }
}

//==============================================================================

    } // namespace Transfer
} // namespace Mutation
