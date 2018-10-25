/*
 * Copyright 2018 von Karman Institute for Fluid Dynamics (VKI)
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

#include "mutation++.h"
#include <catch/catch.hpp>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Kinetics;
using namespace Mutation::Utilities::IO;
using namespace Catch;

void checkReactionType(const std::string& formula, ReactionType type)
{
    // Setup a shared thermo object for testing reaction types
    static Thermodynamics::Thermodynamics thermo(
        "e- N Ar(1) Ar(2) N+ O O- O+ NO NO+ N2 N2+ O2", "RRHO", "ChemNonEq1T");

    XmlElement node(
        "<reaction formula=\"" + formula + "\" >"
        "    <arrhenius A=\"5.0E15\" n=\"+0.0\" T=\"75500.0\" />"
        "</reaction>"
    );

    Reaction reaction(node, thermo);

    INFO("Reaction: " << formula)
    INFO("Should be type: " << reactionTypeString(type))
    INFO("Actual type is: " << reactionTypeString(reaction.type()))

    CHECK(reaction.type() == type);
}

void checkReactionType(
    const std::string& reactants,
    const std::string& products, 
    ReactionType forward_type,
    ReactionType reverse_type)
{
    checkReactionType(reactants + " = " + products, forward_type);
    checkReactionType(products + " = " + reactants, reverse_type);
}

/**
 * Tests if supported reaction types are correctly identified.
 */
TEST_CASE
(
    "Reaction types are correctly identified",
    "[kinetics]"
)
{
    // ASSOCIATIVE_IONIZATION and DISSOCIATIVE_RECOMBINATION
    checkReactionType("N + O", "NO+ + e-", 
        ASSOCIATIVE_IONIZATION, DISSOCIATIVE_RECOMBINATION);
    checkReactionType("N + N", "N2+ + e-", 
        ASSOCIATIVE_IONIZATION, DISSOCIATIVE_RECOMBINATION);

    // CHARGE_EXCHANGE
    checkReactionType("N+ + O = N + O+", CHARGE_EXCHANGE);

    // DISSOCIATION_E and RECOMBINATION_E
    checkReactionType("NO + e-", "N + O + e-", DISSOCIATION_E, RECOMBINATION_E);
    checkReactionType("N2 + e-", "2N + e-",    DISSOCIATION_E, RECOMBINATION_E);

    // DISSOCIATION_M and RECOMBINATION_M
    checkReactionType("NO + M", "N + O + M", DISSOCIATION_M, RECOMBINATION_M);
    checkReactionType("N2 + M", "2N + M",    DISSOCIATION_M, RECOMBINATION_M);
    checkReactionType("NO + N", "2N + O",    DISSOCIATION_M, RECOMBINATION_M);
    checkReactionType("N2 + N", "3N",        DISSOCIATION_M, RECOMBINATION_M);
    
    // ELECTRONIC_ATTACHMENT and ELECTRONIC_DETACHMENT
    // checkReactionType("O + e- + N", "O- + N", 
    //     ELECTRONIC_ATTACHMENT, ELECTRONIC_DETACHMENT);
    // checkReactionType("O + e- + M", "O- + M", 
    //     ELECTRONIC_ATTACHMENT, ELECTRONIC_DETACHMENT);

    // EXCHANGE
    checkReactionType("N2 + O2 = 2NO",   EXCHANGE);
    checkReactionType("NO + O = O2 + N", EXCHANGE);

    // IONIZATION_E and ION_RECOMBINATION_E
    checkReactionType("N + e-", "N+ + 2e-",   IONIZATION_E, ION_RECOMBINATION_E);

    // IONIZATION_M and ION_RECOMBINATION_M
    checkReactionType("N + M", "N+ + M + e-", IONIZATION_M, ION_RECOMBINATION_M);
    checkReactionType("N + N", "N+ + N + e-", IONIZATION_M, ION_RECOMBINATION_M);

    // EXCITATION_E and DEEXCITATION_E
    // checkReactionType("Ar(1) + e-", "Ar(2) + e-", EXCITATION_E, DEEXCITATION_E);

    // EXCITATION_M and DEEXCITATION_M
    // checkReactionType("Ar(1) + M", "Ar(2) + M", EXCITATION_M, DEEXCITATION_M);
}
