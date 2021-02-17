/**
 * @file CollisionPair.h
 *
 * Provides the CollisionPair type.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef TRANSPORT_COLLISION_PAIR_H
#define TRANSPORT_COLLISION_PAIR_H

#include "CollisionIntegral.h"
#include "SharedPtr.h"
#include "XMLite.h"

#include <string>

// Forward declarations
namespace Mutation { namespace Thermodynamics { class Species; }}

namespace Mutation {
    namespace Transport {

/**
 * Enumerates the different type of collision pairs.
 */
enum CollisionType {
    NEUTRAL_NEUTRAL,
    ELECTRON_NEUTRAL,
    ION_NEUTRAL,
    ATTRACTIVE,
    REPULSIVE
};

/**
 * Encapsulates data corresponding to a particular collision pair.
 */
class CollisionPair
{
public:
    /**
     * Loads the collision pair information provided the species in the pair and
     * the root node of the XML collision database.
     */
    CollisionPair(
        const Mutation::Thermodynamics::Species& s1,
        const Mutation::Thermodynamics::Species& s2,
        const Mutation::Utilities::IO::XmlElement* xml);

    // Getter functions
    const Thermodynamics::Species& sp1() const { return *mp_sp1; }
    const Thermodynamics::Species& sp2() const { return *mp_sp2; }
    const std::string& sp1Name() const;
    const std::string& sp2Name() const;
    const std::string& name() const { return m_string; }
    CollisionType type() const { return m_type; }

    /// Get the collision integral corresponding to the given type.
    SharedPtr<CollisionIntegral> get(const std::string& type);

    /// Looks for the XML element representing this pair in the database.
    Mutation::Utilities::IO::XmlElement::const_iterator findPair() const;

private:

    /**
     * Initializes the species names and the type of collision integral this is
     * from the two species objects.
     */
    void initSpeciesData(
        const Mutation::Thermodynamics::Species& s1,
        const Mutation::Thermodynamics::Species& s2);

    /**
     * Returns the iterator pointing to the XmlElement which holds the collision
     * integral of type kind for this collision pair.  If this doesn't exist,
     * then database.end() is returned.
     */
    Mutation::Utilities::IO::XmlElement::const_iterator
    findXmlElementWithIntegralType(const std::string& kind) const;

    /**
     * Loads a particular collision integral from the database.
     */
    SharedPtr<CollisionIntegral> loadIntegral(const std::string& type);

private:

    /// Type of collision
    CollisionType m_type;

    // Species pointers
    const Thermodynamics::Species* mp_sp1;
    const Thermodynamics::Species* mp_sp2;

    // Reference to xml database
    const Mutation::Utilities::IO::XmlElement* mp_xml;

    /// Group of collision integrals loaded for this pair
    std::map<std::string, SharedPtr<CollisionIntegral> > m_integrals;

    /// String representation of this pair
    std::string m_string;
};

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_COLLISION_PAIR_H
