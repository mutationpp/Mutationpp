/**
 * @file CollisionDB.h
 *
 * @brief Declaration of CollisionDB type.
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

#ifndef TRANSPORT_COLLISION_DB_H
#define TRANSPORT_COLLISION_DB_H

#include "CollisionPair.h"
#include "CollisionGroup.h"
#include "XMLite.h"

#include <cassert>
#include <map>
#include <string>
#include <vector>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
    namespace Transport {

/**
 * Provides collision integral data loaded from a database file.
 */
class CollisionDBNew
{
public:

    /**
     * Constructs a new CollisionDB type.
     */
    CollisionDBNew(
        const std::string& db_name,
        const Thermodynamics::Thermodynamics& thermo);

    /// Returns the number of collision pairs in this database.
    int size() const { return m_pairs.size(); }

    /// Returns i'th CollisionPair in the database.
    const CollisionPairNew& operator [] (int i) const {
        assert(i >= 0 && i < size());
        return m_pairs[i];
    }

    /**
     * Returns the collision group corresponding to name, updated to the current
     * state.  If the group is not loaded yet, it is first loaded from the
     * database.
     */
    const CollisionGroup& group(const std::string& name);

    /// Provides Q11 collision integrals for electron-heavy interactions.
    const CollisionGroup& Q11ei() { return group("Q11ei"); }

    /// Provides Q11 collision integrals for heavy-heavy interactions.
    const CollisionGroup& Q11ij() { return group("Q11ij"); }

    /// Provides Q22 collision integrals for electron-heavy interactions.
    const CollisionGroup& Q22ei() { return group("Q22ei"); }

    /// Provides Q22 collision integrals for heavy-heavy interactions.
    const CollisionGroup& Q22ij() { return group("Q22ij"); }

    /// Provides B* collision integrals for electron-heavy interactions.
    const CollisionGroup& Bstei() { return group("Bstei"); }

    /// Provides B* collision integrals for heavy-heavy interactions.
    const CollisionGroup& Bstij() { return group("Bstij"); }

    /// Provides C* collision integrals for electron-heavy interactions.
    const CollisionGroup& Cstei() { return group("Cstei"); }

    /// Provides C* collision integrals for heavy-heavy interactions.
    const CollisionGroup& Cstij() { return group("Cstij"); }

private:

    Mutation::Utilities::IO::XmlDocument m_database;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    bool m_tabulate;

    // List of collision pairs
    std::vector<CollisionPairNew> m_pairs;

    // CollisionGroup container
    std::map<std::string, CollisionGroup> m_groups;
};

    } // namespace Transport
} // namespace Mutation


#endif // TRANSPORT_COLLISION_DB_H
