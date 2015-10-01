/**
 * @file CollisionDB.cpp
 *
 * @brief Implementation of CollisionDB type.
 */

/*
 * Copyright 2014-2015 von Karman Institute for Fluid Dynamics (VKI)
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

#include "CollisionDBNew.h"
#include "Species.h"
#include "Thermodynamics.h"
#include "XMLite.h"
#include "Utilities.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities::IO;

#include <fstream>
using namespace std;

namespace Mutation {
    namespace Transport {

//==============================================================================

CollisionDBNew::CollisionDBNew(
    const string& db_name, const Thermodynamics::Thermodynamics& thermo) :
    m_thermo(thermo)
{
    // First load all of the necessary collision pairs
    loadCollisionPairs(db_name);

    // Fill the CollisionGroups
    const int k = (thermo.hasElectrons() ? thermo.nSpecies() : 0);
    // Electron-Heavy
    vector<CollisionPairNew>::const_iterator start = m_pairs.begin();
    vector<CollisionPairNew>::const_iterator end   = start+k;
    m_Q11ei.manage(start, end, &CollisionPairNew::Q11);
    m_Q22ei.manage(start, end, &CollisionPairNew::Q22);
    m_Bstei.manage(start, end, &CollisionPairNew::Bst);
    m_Cstei.manage(start, end, &CollisionPairNew::Cst);
    // Heavy-Heavy
    start = end;
    end   = m_pairs.end();
    m_Q11ij.manage(start, end, &CollisionPairNew::Q11);
    m_Q22ij.manage(start, end, &CollisionPairNew::Q22);
    m_Bstij.manage(start, end, &CollisionPairNew::Bst);
    m_Cstij.manage(start, end, &CollisionPairNew::Cst);
}

//==============================================================================

void CollisionDBNew::loadCollisionPairs(string db_name)
{
    // Add extension if necessary
    if (db_name.substr(db_name.length()-5) != ".xml")
        db_name = db_name + ".xml";

    // Add directory
    {
        ifstream file(db_name.c_str(), ios::in);

        // If that doesn't work, look in MPP_DATA_DIRECTORY/transport
        if (!file.is_open())
            db_name = Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") +
                "/transport/" + db_name;
    }

    // Now load the XML file
    XmlDocument db_doc(db_name);
    XmlElement root = db_doc.root();

    // Loop over the species and create the list of species pairs
    const vector<Species>& species = m_thermo.species();
    for (int i = 0; i < species.size(); ++i)
        for (int j = i; j < species.size(); ++j)
            m_pairs.push_back(CollisionPairNew(species[i], species[j], root));
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Q11ei() {
    return m_Q11ei.update(m_thermo.Te(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Q11ij() {
    return m_Q11ij.update(m_thermo.T(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Q22ei() {
    return m_Q22ei.update(m_thermo.Te(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Q22ij() {
    return m_Q22ij.update(m_thermo.T(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Bstei() {
    return m_Bstei.update(m_thermo.Te(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Bstij() {
    return m_Bstij.update(m_thermo.T(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Cstei() {
    return m_Cstei.update(m_thermo.Te(), m_thermo);
}

//==============================================================================

const CollisionGroup& CollisionDBNew::Cstij() {
    return m_Cstij.update(m_thermo.T(), m_thermo);
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation

