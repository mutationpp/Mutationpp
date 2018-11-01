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

#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include <catch.hpp>
#include <Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Catch;
using namespace Eigen;


/**
 * Tests the SpeciesNameFSM::parse() method.
 */
TEST_CASE("Parsing species name with SpeciesNameFSM", "[thermodynamics]")
{
    SpeciesNameFSM sm;

    CHECK(sm.parse("e-"));
    CHECK(sm.stoichiometry().size() == 1);
    CHECK(sm.stoichiometry().find("e-")->second == 1);

    CHECK(sm.parse("ABC"));
    CHECK(sm.stoichiometry().size() == 3);
    CHECK(sm.stoichiometry().find("A")->second == 1);
    CHECK(sm.stoichiometry().find("B")->second == 1);
    CHECK(sm.stoichiometry().find("C")->second == 1);

    CHECK(sm.parse("ABC+"));
    CHECK(sm.stoichiometry().size() == 4);
    CHECK(sm.stoichiometry().find("A")->second == 1);
    CHECK(sm.stoichiometry().find("B")->second == 1);
    CHECK(sm.stoichiometry().find("C")->second == 1);
    CHECK(sm.stoichiometry().find("e-")->second == -1);

    CHECK(sm.parse("A-"));
    CHECK(sm.stoichiometry().size() == 2);
    CHECK(sm.stoichiometry().find("A")->second == 1);
    CHECK(sm.stoichiometry().find("e-")->second == 1);

    CHECK(sm.parse("AaBbCc"));
    CHECK(sm.stoichiometry().size() == 3);
    CHECK(sm.stoichiometry().find("Aa")->second == 1);
    CHECK(sm.stoichiometry().find("Bb")->second == 1);
    CHECK(sm.stoichiometry().find("Cc")->second == 1);

    CHECK(sm.parse("Aa1B3Cc1++"));
    CHECK(sm.stoichiometry().size() == 4);
    CHECK(sm.stoichiometry().find("e-")->second == -2);
    CHECK(sm.stoichiometry().find("Aa")->second == 1);
    CHECK(sm.stoichiometry().find("B")->second == 3);
    CHECK(sm.stoichiometry().find("Cc")->second == 1);

    CHECK(sm.parse("C6H5OH+"));
    CHECK(sm.stoichiometry().size() == 4);
    CHECK(sm.stoichiometry().find("e-")->second == -1);
    CHECK(sm.stoichiometry().find("C")->second == 6);
    CHECK(sm.stoichiometry().find("H")->second == 6);
    CHECK(sm.stoichiometry().find("O")->second == 1);

    CHECK_FALSE(sm.parse(""));
    CHECK_FALSE(sm.parse(" "));
    CHECK_FALSE(sm.parse("23A"));
    CHECK_FALSE(sm.parse("!e#5Hz"));
    CHECK_FALSE(sm.parse("Aaa"));
    CHECK_FALSE(sm.parse("A+B"));

    CHECK(sm.stoichiometry().size() == 0);
}


/**
 * Tests the initialization of Species objects using different constructors and
 * the assignment operator.
 */
TEST_CASE("Initialization of Species objects", "[thermodynamics]")
{
    const double mw_N = 0.0140067;
    const double mw_O = 0.0159994;
    const double mw_e = 0.00000055;

    // Initialize from species name
    CHECK_NOTHROW(Species("NO2+"));
    Species s("NO2+");
    CHECK(s.charge() == 1);
    CHECK(s.groundStateName() == "NO2+");
    CHECK(s.isIon());
    CHECK(s.level() == 0);
    CHECK(s.molecularWeight() == mw_N + 2*mw_O - mw_e);
    CHECK(s.nAtoms("N") == 1);
    CHECK(s.nAtoms("O") == 2);
    CHECK(s.nAtoms("e-") == -1);
    CHECK(s.nAtoms() == 3);
    CHECK(s.name() == "NO2+");
    CHECK(s.phase() == GAS);
    CHECK(s.stoichiometry()[0] == std::make_pair(std::string("N"), 1));
    CHECK(s.stoichiometry()[1] == std::make_pair(std::string("O"), 2));
    CHECK(s.stoichiometry()[2] == std::make_pair(std::string("e-"), -1));
    CHECK(s.type() == MOLECULE);

    // Assignment operator and initialize from name + phase
    CHECK_NOTHROW(s = Species("O", SOLID));
    CHECK(s.charge() == 0);
    CHECK(s.groundStateName() == "O");
    CHECK_FALSE(s.isIon());
    CHECK(s.level() == 0);
    CHECK(s.molecularWeight() == mw_O);
    CHECK(s.nAtoms("O") == 1);
    CHECK(s.nAtoms() == 1);
    CHECK(s.name() == "O");
    CHECK(s.phase() == SOLID);
    CHECK(s.stoichiometry()[0] == std::make_pair(std::string("O"), 1));
    CHECK(s.type() == ATOM);

    // Assignment and initialize from name + phase + stoichiometry
    CHECK_NOTHROW(s = Species(
         "Nitrogen Ion", GAS, Species::StoichList()("N", 1)("e-", -1)));
    CHECK(s.charge() == 1);
    CHECK(s.groundStateName() == "Nitrogen Ion");
    CHECK(s.isIon());
    CHECK(s.level() == 0);
    CHECK(s.molecularWeight() == mw_N - mw_e);
    CHECK(s.nAtoms("N") == 1);
    CHECK(s.nAtoms("e-") == -1);
    CHECK(s.nAtoms() == 1);
    CHECK(s.name() == "Nitrogen Ion");
    CHECK(s.phase() == GAS);
    CHECK(s.stoichiometry()[0] == std::make_pair(std::string("N"), 1));
    CHECK(s.stoichiometry()[1] == std::make_pair(std::string("e-"), -1));
    CHECK(s.type() == ATOM);

    // Copy constructor
    s = Species(s);
    CHECK(s.charge() == 1);
    CHECK(s.groundStateName() == "Nitrogen Ion");
    CHECK(s.isIon());
    CHECK(s.level() == 0);
    CHECK(s.molecularWeight() == mw_N - mw_e);
    CHECK(s.nAtoms("N") == 1);
    CHECK(s.nAtoms("e-") == -1);
    CHECK(s.nAtoms() == 1);
    CHECK(s.name() == "Nitrogen Ion");
    CHECK(s.phase() == GAS);
    CHECK(s.stoichiometry()[0] == std::make_pair(std::string("N"), 1));
    CHECK(s.stoichiometry()[1] == std::make_pair(std::string("e-"), -1));
    CHECK(s.type() == ATOM);

    // Initialize level from species
    s = Species(s, 3);
    CHECK(s.charge() == 1);
    CHECK(s.groundStateName() == "Nitrogen Ion");
    CHECK(s.isIon());
    CHECK(s.level() == 3);
    CHECK(s.molecularWeight() == mw_N - mw_e);
    CHECK(s.nAtoms("N") == 1);
    CHECK(s.nAtoms("e-") == -1);
    CHECK(s.nAtoms() == 1);
    CHECK(s.name() == "Nitrogen Ion(3)");
    CHECK(s.phase() == GAS);
    CHECK(s.stoichiometry()[0] == std::make_pair(std::string("N"), 1));
    CHECK(s.stoichiometry()[1] == std::make_pair(std::string("e-"), -1));
    CHECK(s.type() == ATOM);

    std::stringstream ss;
    ss << "<species name=\"N2\">\n";
    ss << "    <stoichiometry> N : 2 </stoichiometry>\n";
    ss << "</species>";

    // Initialize from XML
    CHECK_NOTHROW(s = Species(Utilities::IO::XmlElement(ss.str())));
    CHECK(s.charge() == 0);
    CHECK(s.groundStateName() == "N2");
    CHECK_FALSE(s.isIon());
    CHECK(s.level() == 0);
    CHECK(s.molecularWeight() == 2*mw_N);
    CHECK(s.nAtoms("N") == 2);
    CHECK(s.nAtoms() == 2);
    CHECK(s.name() == "N2");
    CHECK(s.phase() == GAS);
    CHECK(s.stoichiometry()[0] == std::make_pair(std::string("N"), 2));
    CHECK(s.type() == MOLECULE);

    CHECK_THROWS_AS(Species("N+O"), InvalidInputError);
    CHECK_THROWS_AS(Species("AbN2"), InvalidInputError);
}



