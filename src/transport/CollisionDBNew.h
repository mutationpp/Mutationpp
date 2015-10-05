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

    /// Returns number of species in the database.
    int nSpecies() const;

    /// Returns number of heavy particles in the database.
    int nHeavy() const;

    /// Returns the number of collision pairs in this database.
    int size() const { return m_pairs.size(); }

    /// Returns an array of species masses in kg.
    const Eigen::ArrayXd& mass() const { return m_mass; }

    /// Returns the reference of thermo.
    const Thermodynamics::Thermodynamics& thermo() const { return m_thermo; }

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
    const Eigen::ArrayXd& Q11ei() { return group("Q11ei").array(); }

    /// Provides Q11 collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Q11ij() { return group("Q11ij").array(); }

    /// Provides Q12 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q12ei() { return group("Q12ei").array(); }

    /// Provides Q13 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q13ei() { return group("Q13ei").array(); }

    /// Provides Q14 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q14ei() { return group("Q14ei").array(); }

    /// Provides Q15 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q15ei() { return group("Q15ei").array(); }

    /// Provides Q22 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q22ei() { return group("Q22ei").array(); }

    /// Provides Q22 collision integrals for diagonal heavy-heavy interactions.
    const Eigen::ArrayXd& Q22ii() { return group("Q22ii").array(); }

    /// Provides Q22 collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Q22ij() { return group("Q22ij").array(); }

    /// Provides Q23 collision integrals for electron-electron interaction.
    double Q23ee() { return group("Q23ee")[0]; }

    /// Provides Q24 collision integrals for electron-electron interaction.
    double Q24ee() { return group("Q24ee")[0]; }

    /// Provides A* collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Astei() { return group("Astei").array(); }

    /// Provides A* collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Astij() { return group("Astij").array(); }

    /// Provides B* collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Bstei() { return group("Bstei").array(); }

    /// Provides B* collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Bstij() { return group("Bstij").array(); }

    /// Provides C* collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Cstei() { return group("Cstei").array(); }

    /// Provides C* collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Cstij() { return group("Cstij").array(); }

    /// Returns the pure, heavy species shear viscosities array.
    const Eigen::ArrayXd& etai();

    /// Returns binary diffusion coefficients for electron-heavy interactions.
    const Eigen::ArrayXd& nDei();

    /// Returns binary diffusion coefficients for heavy-heavy interactions.
    const Eigen::ArrayXd& nDij();

    /**
     * Returns the average diffusion coefficients.
     * \f[ D_{im} = \frac{(1-x_i)}{\sum_{j\ne i}x_j/\mathscr{D}_{ij}} \f]
     */
    const Eigen::ArrayXd& Dim();

private:

    /// Types of collision integral groups (ie: electron-heavy, ...).
    enum GroupType {
        EE = 0, EI, II, IJ, BAD_TYPE
    };

    /// Determines the type of group from the group name.
    GroupType groupType(const std::string& name);

private:

    Mutation::Utilities::IO::XmlDocument m_database;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    // Tabulation parameters
    bool m_tabulate;
    double m_table_min;
    double m_table_max;
    double m_table_del;

    // List of collision pairs
    std::vector<CollisionPairNew> m_pairs;

    // CollisionGroup container
    std::map<std::string, CollisionGroup> m_groups;

    // Derived data and helper arrays
    Eigen::ArrayXd m_mass;
    Eigen::ArrayXd m_etai;
    Eigen::ArrayXd m_etafac;
    Eigen::ArrayXd m_nDei;
    Eigen::ArrayXd m_Deifac;
    Eigen::ArrayXd m_nDij;
    Eigen::ArrayXd m_Dijfac;
    Eigen::ArrayXd m_Dim;
};

    } // namespace Transport
} // namespace Mutation


#endif // TRANSPORT_COLLISION_DB_H
