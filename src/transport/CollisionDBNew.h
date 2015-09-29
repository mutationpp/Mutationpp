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

#include "Species.h"
#include "CollisionPair.h"

#include <cassert>
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

private:

    void loadCollisionPairs(
        std::string db_name,
        const std::vector<Thermodynamics::Species>& species);

private:

    std::vector<CollisionPairNew> m_pairs;
};

    } // namespace Transport
} // namespace Mutation


#endif // TRANSPORT_COLLISION_DB_H
