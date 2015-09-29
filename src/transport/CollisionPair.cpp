/**
 * @file CollisionPair.cpp
 *
 * Implements the CollisionPair type.
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

#include "AutoRegistration.h"
#include "CollisionPair.h"
#include "Species.h"
#include "XMLite.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

#include <string>
#include <iostream> // remove
using namespace std;

namespace Mutation {
    namespace Transport {

//==============================================================================

CollisionPairNew::CollisionPairNew(
    const Species& s1, const Species& s2, const IO::XmlElement& xml)
{
    // First initialize species info
    initSpeciesData(s1, s2);

    // Next load up the collision integrals
    mp_Q11 = loadIntegral("Q11", xml);
    mp_Q22 = loadIntegral("Q22", xml);
    mp_Bst = loadIntegral("Bst", xml);
    mp_Cst = loadIntegral("Cst", xml);
}

//==============================================================================

void CollisionPairNew::initSpeciesData(const Species& s1, const Species& s2)
{
    string name1 = s1.groundStateName();
    string name2 = s2.groundStateName();

    // Determine the type of collision
    if (!s1.isIon()) { // first species is neutral
        if (s2.type() == ELECTRON)
            m_type = ELECTRON_NEUTRAL;
        else if (s2.isIon())
            m_type = ION_NEUTRAL;
        else
            m_type = NEUTRAL_NEUTRAL;
    } else { // first species is an ion (electron or heavy)
        if (s2.isIon()) {
            if (s1.charge() * s2.charge() > 0)
                m_type = REPULSIVE;
            else
                m_type = ATTRACTIVE;
        } else if (s1.type() == ELECTRON)
            m_type = ELECTRON_NEUTRAL;
        else
            m_type = ION_NEUTRAL;
    }

    // Initialize the species names (e- first, then alphabetical order)
    m_sp1 = s1.groundStateName();
    m_sp2 = s2.groundStateName();

    if (m_sp1 > m_sp2)
        std::swap(m_sp1, m_sp2);
    if (m_sp2 == "e-")
        std::swap(m_sp1, m_sp2);
}

//==============================================================================

IO::XmlElement::const_iterator CollisionPairNew::findXmlElementWithIntegralType(
    const string& kind, const IO::XmlElement& database) const
{
    // First check if this collision pair is explicitly given in the database
    IO::XmlElement::const_iterator iter = database.findTag("pair");
    string sp1, sp2;
    while (iter != database.end()) {
        iter->getAttribute("s1", sp1, "Collision pair missing sp1 attribute.");
        iter->getAttribute("s2", sp2, "Collision pair missing sp2 attribute.");

        // Check if species names match
        if ((m_sp1 == sp1 && m_sp2 == sp2) || (m_sp1 == sp2 && m_sp2 == sp1))
            break;

        // Get next collision pair in the database
        iter = database.findTag("pair", ++iter);
    }

    // Found the pair, so check if the integral is explicitly given
    if (iter != database.end()) {
        IO::XmlElement::const_iterator Qiter = iter->findTag(kind);
        if (Qiter != iter->end())
            return Qiter;
    }

    // Didn't find the pair or integral wasn't given, so look in defaults
    iter = database.findTag("defaults");
    if (iter == database.end()) return iter;

    IO::XmlElement::const_iterator pair;
    switch (m_type) {
        case NEUTRAL_NEUTRAL:
            pair = iter->findTag("neutral-neutral"); break;
        case ELECTRON_NEUTRAL:
            pair = iter->findTag("electron-neutral"); break;
        case ION_NEUTRAL:
            pair = iter->findTag("ion-neutral"); break;
        default:
            pair = iter->findTag("charged"); break;
    }
    if (pair == iter->end())
        return database.end();

    iter = pair->findTag(kind);
    if (iter != pair->end())
        return iter;

    return database.end();
}


//==============================================================================

SharedPtr<CollisionIntegral> CollisionPairNew::loadIntegral(
    const string& kind, const IO::XmlElement& database) const
{
    IO::XmlElement::const_iterator iter =
        findXmlElementWithIntegralType(kind, database);
    if (iter == database.end())
        return SharedPtr<CollisionIntegral>();

    string type;
    iter->getAttribute("type", type, "Integral type must be specified.");

    CollisionIntegral* p_ci(
        Config::Factory<CollisionIntegral>::create(
            type, CollisionIntegral::ARGS(*iter, *this)));

    if (p_ci == NULL) iter->parseError(
        "Invalid collision integral type '" + type + "' for " + kind +
		"_(" + m_sp1 + ", " + m_sp2 + ").");

    return SharedPtr<CollisionIntegral>(p_ci);
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
