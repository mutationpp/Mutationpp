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

#ifndef SPECIES_NAME_FSM_H
#define SPECIES_NAME_FSM_H

#include <string>
#include <map>

namespace Mutation {
namespace Thermodynamics {

/**
 * Implements the finite-state machine (FSM) logic for parsing a species name
 * and retrieving the species stoichiometry.
 */
class SpeciesNameFSM
{
public:
    /// Constructor does nothing.
    SpeciesNameFSM();

    /**
     * Parses the given species name and extracts information.  Returns true if
     * the name was properly formatted, false otherwise.
     */
    bool parse(const std::string& name);

    /**
     * Returns the parsed stoichiometry in terms of a map of element strings to
     * integer numbers.
     */
    const std::map<std::string, int>& stoichiometry() const {
        return m_stoich;
    }

private:

    /// Run single character through the FSM.
    void next(char c);

    /// Behavior for state FirstChar.
    void stateFirstChar(char c);

    /// Behavior for state SecondChar.
    void stateSecondChar(char c);

    /// Behavior for state Number.
    void stateNumber(char c);

    /// Behavior for state Charge
    void stateCharge(char c);

    /// Adds n atoms of given element to the stoichiometry map.
    void incrementElement(const std::string& element, int n);

private:

    /// The current state.
    enum {
        FirstChar, SecondChar, Number, Charge, BadFormat
    } m_state;

    /// The current element string
    std::string m_element;

    /// The current number of elements
    int         m_number;

    /// The stoichiometry map
    std::map<std::string, int> m_stoich;

}; // class SpeciesNameFSM

} // namespace Thermodynamics
} // namespace Mutation

#endif // SPECIES_NAME_FSM_H
