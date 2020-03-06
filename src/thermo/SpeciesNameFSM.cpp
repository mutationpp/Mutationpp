/*
 * Copyright 2017-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "SpeciesNameFSM.h"
#include <algorithm>
#include <functional>

namespace Mutation {
namespace Thermodynamics {

//==============================================================================

SpeciesNameFSM::SpeciesNameFSM()
{
    m_state = FirstChar;
    m_element = "";
    m_number = 0;
    m_stoich.clear();
}

//==============================================================================

bool SpeciesNameFSM::parse(const std::string& s)
{
    // Initialize
    m_state = FirstChar;
    m_stoich.clear();

    // Empty string is considered a bad format
    if (s.empty())
        return false;

    // Special case: free electron
    if (s == "e-") {
        m_stoich["e-"] = 1;
        return true;
    }

    // Run each character through the FSM
    for_each(s.begin(), s.end(), std::bind1st(std::mem_fun(&SpeciesNameFSM::next), this));
    next('\0');

    if (m_state == BadFormat) {
        m_stoich.clear();
        return false;
    }

    return true;
}

//==============================================================================

void SpeciesNameFSM::next(char c)
{
    switch(m_state) {
    case FirstChar:  stateFirstChar(c);  return;
    case SecondChar: stateSecondChar(c); return;
    case Number:     stateNumber(c);     return;
    case Charge:     stateCharge(c);     return;
    case BadFormat:                      return;
    }
}

//==============================================================================

void SpeciesNameFSM::stateFirstChar(char c)
{
    // Check end state
    if (c == '\0')
        return;

    // Expecting [A-Z]
    if (c >= 'A' && c <= 'Z') {
        m_element = c;
        m_state = SecondChar;
    // Anything else is an error
    } else
        m_state = BadFormat;
}

//==============================================================================

void SpeciesNameFSM::stateSecondChar(char c)
{
    // Check end state
    if (c == '\0') {
        incrementElement(m_element, 1);
        return;
    }

    // Expecting [a-z]
    if (c >= 'a' && c <= 'z') {
        m_element += c;
        m_number = 0;
        m_state = Number;
    // If this is actually a number then go to state Number
    } else if (c >= '0' && c <= '9') {
        m_state = Number;
        m_number = 0;
        stateNumber(c);
    // A capital letter indicates we are starting a new element
    } else if (c >= 'A' && c <= 'Z') {
        m_state = FirstChar;
        incrementElement(m_element, 1);
        stateFirstChar(c);
    // A plus or minus indicates charge
    } else if (c == '+' || c == '-') {
        m_state = Charge;
        incrementElement(m_element, 1);
        stateCharge(c);
    // Otherwise we have an error
    } else
        m_state = BadFormat;
}

//==============================================================================

void SpeciesNameFSM::stateNumber(char c)
{
    // Check end state
    if (c == '\0') {
        incrementElement(m_element, std::max(m_number, 1));
        return;
    }

    // Expecting a digit
    if (c >= '0' && c <= '9') {
        m_number = 10*m_number + (c - '0');
    // If we have a capital letter, then we are starting a new element
    } else if (c >= 'A' && c <= 'Z') {
        m_state = FirstChar;
        incrementElement(m_element, std::max(m_number, 1));
        stateFirstChar(c);
    // A plus or minus indicates charge
    } else if (c == '+' || c == '-') {
        m_state = Charge;
        incrementElement(m_element, std::max(m_number, 1));
        stateCharge(c);
    // Otherwise we have an error
    } else
        m_state = BadFormat;
}

//==============================================================================

void SpeciesNameFSM::stateCharge(char c)
{
    // Check end state
    if (c == '\0')
        return;

    // '+' decreases electron
    if (c == '+')
        incrementElement("e-", -1);
    // '-' increases electron
    else if (c == '-')
        incrementElement("e-", 1);
    // Anything else is considered an error
    else
        m_state = BadFormat;
}

//==============================================================================

void SpeciesNameFSM::incrementElement(const std::string& element, int n)
{
    if (m_stoich.find(element) == m_stoich.end())
        m_stoich[element] = n;
    else
        m_stoich[element] += n;
}

//==============================================================================

} // namespace Thermodynamics
} // namespace Mutation
