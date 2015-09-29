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
    const string& db_name, const Thermodynamics::Thermodynamics& thermo)
{
    loadCollisionPairs(db_name, thermo.species());
}

//==============================================================================

void CollisionDBNew::loadCollisionPairs(
    string db_name, const vector<Species>& species)
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
    for (int i = 0; i < species.size(); ++i)
        for (int j = i; j < species.size(); ++j)
            m_pairs.push_back(CollisionPairNew(species[i], species[j], root));
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation

