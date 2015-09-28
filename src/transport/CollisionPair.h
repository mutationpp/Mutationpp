/**
 * @file CollisionPair.h
 *
 * Provides the CollisionPair type.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015 von Karman Institute for Fluid Dynamics (VKI)
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

#include <string>

namespace Mutation {
    namespace Thermodynamics { class Species; }
    namespace Utilities {
        namespace IO { class XmlElement; }
    }

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

// Forward declarations
class CollisionIntegral;

/**
 * Encapsulates data corresponding to a particular collision pair.
 */
class CollisionPairNew
{
public:
    /**
     * Loads the collision pair information provided the species in the pair and
     * the root node of the XML collision database.
     */
    CollisionPairNew(
        const Thermodynamics::Species& s1, const Thermodynamics::Species& s2,
        const Utilities::IO::XmlElement& xml);

    /**
     * Destructor.
     */
    ~CollisionPairNew();

    // Getter functions
    const std::string& species1() const { return m_sp1; }
    const std::string& species2() const { return m_sp2; }

    const CollisionIntegral* const Q11() { return mp_Q11; }
    const CollisionIntegral* const Q22() { return mp_Q22; }
    const CollisionIntegral* const Bst() { return mp_Bst; }
    const CollisionIntegral* const Cst() { return mp_Cst; }

private:

    /**
     * Initializes the species names and the type of collision integral this is
     * from the two species objects.
     */
    void initSpeciesData(
        const Thermodynamics::Species& s1,  const Thermodynamics::Species& s2);

    /**
     * Returns the iterator pointing to the XmlElement which holds the collision
     * integral of type kind for this collision pair.  If this doesn't exist,
     * then database.end() is returned.
     */
    Utilities::IO::XmlElement::const_iterator
    findXmlElementWithIntegralType(
        const string& kind, const XmlElement& database);

    /**
     * Loads a particular collision integral from the database.
     */
    CollisionIntegral* loadIntegral(
        const std::string& type, const Utilities::IO::XmlElement& database);

private:

    CollisionType m_type;
    std::string   m_sp1;
    std::string   m_sp2;

    // collision integrals
    CollisionIntegral* mp_Q11;
    CollisionIntegral* mp_Q22;
    CollisionIntegral* mp_Bst;
    CollisionIntegral* mp_Cst;
};

    } // namespace Transport
} // namespace Mutation
