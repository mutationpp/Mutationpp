/**
 * @file CollisionPair.cpp
 *
 * Implements the CollisionPair type.
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

#include "AutoRegistration.h"
#include "CollisionPair.h"
#include "Species.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

#include <string>
#include <iostream> // remove
using namespace std;

namespace Mutation {
    namespace Transport {

//==============================================================================

CollisionPair::CollisionPair(
    const Species& s1, const Species& s2, const IO::XmlElement* xml) :
    mp_xml(xml)
{
    // First initialize species info
    initSpeciesData(s1, s2);
}

//==============================================================================

const string& CollisionPair::sp1Name() const {
    return mp_sp1->groundStateName();
}

//==============================================================================

const string& CollisionPair::sp2Name() const {
    return mp_sp2->groundStateName();
}

//==============================================================================

SharedPtr<CollisionIntegral> CollisionPair::get(const string& type)
{
    map<string, SharedPtr<CollisionIntegral> >::iterator iter =
        m_integrals.find(type);
    if (iter != m_integrals.end())
        return iter->second;

    return m_integrals.insert(
        make_pair(type, loadIntegral(type))).first->second;
}

//==============================================================================

void CollisionPair::initSpeciesData(const Species& s1, const Species& s2)
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

    // Store species in order depending on their names (e- always first)
    mp_sp1 = &s1;
    mp_sp2 = &s2;

    if (sp1Name() > sp2Name())
        std::swap(mp_sp1, mp_sp2);
    if (sp2Name() == "e-")
        std::swap(mp_sp1, mp_sp2);

    // Make the string representation of this pair
    m_string = "(" + sp1Name() + "," + sp2Name() + ")";
}

//==============================================================================

IO::XmlElement::const_iterator CollisionPair::findPair() const
{
    IO::XmlElement::const_iterator iter = mp_xml->findTag("pair");
    string sp1, sp2;
    while (iter != mp_xml->end()) {
        iter->getAttribute("s1", sp1, "Collision pair missing sp1 attribute.");
        iter->getAttribute("s2", sp2, "Collision pair missing sp2 attribute.");

        // Check if species names match
        if ((sp1Name() == sp1 && sp2Name() == sp2) ||
            (sp1Name() == sp2 && sp2Name() == sp1))
            break;

        // Get next collision pair in the database
        iter = mp_xml->findTag("pair", ++iter);
    }

    return iter;
}

//==============================================================================

IO::XmlElement::const_iterator
CollisionPair::findXmlElementWithIntegralType(
    const string& kind) const
{
    // First check if this collision pair is explicitly given in the database
    IO::XmlElement::const_iterator iter = findPair();

    // Found the pair, so check if the integral is explicitly given
    if (iter != mp_xml->end()) {
        IO::XmlElement::const_iterator Qiter = iter->findTag(kind);
        if (Qiter != iter->end())
            return Qiter;
    }

    // Didn't find the pair or integral wasn't given, so look in defaults
    iter = mp_xml->findTag("defaults");
    if (iter == mp_xml->end()) return iter;

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
        return mp_xml->end();

    iter = pair->findTag(kind);
    if (iter != pair->end())
        return iter;

    return mp_xml->end();
}


//==============================================================================

SharedPtr<CollisionIntegral> CollisionPair::loadIntegral(const string& kind)
{
    IO::XmlElement::const_iterator iter =
        findXmlElementWithIntegralType(kind);

    if (iter == mp_xml->end()) {
        throw Mutation::MissingDataError() 
            << "Collision integral " << kind << " is not given for the pair "
            << name() << '.';
    }

    return CollisionIntegral::load(CollisionIntegral::ARGS(*iter, *this, kind));
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
