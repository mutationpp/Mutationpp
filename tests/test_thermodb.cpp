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
 * Checks that species list descriptors are correctly parsed.
 */
TEST_CASE("SpeciesListDescriptor tests", "[thermodynamics][loading]")
{
    // Hard coded list of available species
    std::vector<Species> available;
    available.push_back(Species("N"));             // 0: N
    available.push_back(Species(available[0], 0)); // 1: N(0)
    available.push_back(Species(available[0], 1)); // 2: N(1)
    available.push_back(Species(available[0], 2)); // 3: N(2)
    available.push_back(Species("O"));             // 4: O
    available.push_back(Species("O2"));            // 5: O2
    available.push_back(Species("NO"));            // 6: NO
    available.push_back(Species("N2"));            // 7: N2
    available.push_back(Species("Ar"));            // 8: Ar
    available.push_back(                           // 9: C(gr)
        Species("C(gr)", SOLID, Species::StoichList()("C", 1)));
    available.push_back(Species("e-"));            // 10: e-

    SECTION("Check matches and ordering without implicit rules") {
        SpeciesListDescriptor ld("C(gr) N O e- \"N2\" NO O2");
        CHECK(ld.matches(available[0])); // N
        CHECK_FALSE(ld.matches(available[1])); // N(0)
        CHECK_FALSE(ld.matches(available[2])); // N(1)
        CHECK_FALSE(ld.matches(available[3])); // N(2)
        CHECK(ld.matches(available[4])); // O
        CHECK(ld.matches(available[5])); // O2
        CHECK(ld.matches(available[6])); // NO
        CHECK(ld.matches(available[7])); // N2
        CHECK_FALSE(ld.matches(available[8])); // Ar
        CHECK(ld.matches(available[9])); // C(gr)
        CHECK(ld.matches(available[10])); // e-

        std::list<Species> input;
        std::vector<Species>::const_iterator it;
        for (it = available.begin(); it != available.end(); ++it)
            if (ld.matches(*it)) input.push_back(*it);

        std::vector<Species> output;
        std::vector<std::string> missing;
        ld.order(input, output, missing);

        // Correct ordering: e- N O N2 NO O2 C(gr)
        CHECK(output.size() == 7);
        CHECK(output[0].name() == "e-");
        CHECK(output[1].name() == "N");
        CHECK(output[2].name() == "O");
        CHECK(output[3].name() == "N2");
        CHECK(output[4].name() == "NO");
        CHECK(output[5].name() == "O2");
        CHECK(output[6].name() == "C(gr)");
        CHECK(missing.size() == 0);
    }

    SECTION("Check implicit rules") {
        SpeciesListDescriptor ld("N(*) {solids with C} \"O\"");
        CHECK_FALSE(ld.matches(available[0])); // N
        CHECK(ld.matches(available[1])); // N(0)
        CHECK(ld.matches(available[2])); // N(1)
        CHECK(ld.matches(available[3])); // N(2)
        CHECK(ld.matches(available[4])); // O
        CHECK_FALSE(ld.matches(available[5])); // O2
        CHECK_FALSE(ld.matches(available[6])); // NO
        CHECK_FALSE(ld.matches(available[7])); // N2
        CHECK_FALSE(ld.matches(available[8])); // Ar
        CHECK(ld.matches(available[9])); // C(gr)
        CHECK_FALSE(ld.matches(available[10])); // e-

        std::list<Species> input;
        std::vector<Species>::const_iterator it;
        for (it = available.begin(); it != available.end(); ++it)
            if (ld.matches(*it)) input.push_back(*it);

        std::vector<Species> output;
        std::vector<std::string> missing;
        ld.order(input, output, missing);

        // Correct ordering: N(0) N(1) N(2) O C(gr)
        CHECK(output.size() == 5);
        CHECK(output[0].name() == "N(0)");
        CHECK(output[1].name() == "N(1)");
        CHECK(output[2].name() == "N(2)");
        CHECK(output[3].name() == "O");
        CHECK(output[4].name() == "C(gr)");
        CHECK(missing.size() == 0);
    }
}

// Flags for checkThermoDBLoad()
int DEFAULTS   = 1;
int AIR5       = 2;
int MULTIPHASE = 4;

/**
 * Performs generic tests on a ThermoDB* type to see if it can properly load
 * species thermodynamic data.  It is assumed that the species N O NO N2 O2 are
 * present in the database.  For multiphase checks, C(gr) C CO CO2 should also
 * be present.
 *
 * @param db The ThermoDB* object
 * @param check_multiphase true if multiphase should be checked
 */
void checkThermoDBLoad(ThermoDB* db, int flags)
{
    const double mw_N = 0.0140067;
    const double mw_O = 0.0159994;
    const double mw_C = 0.012011;

    if (flags & DEFAULTS) {
        SECTION("defaults") {
            CHECK(db->standardPressure() == ONEATM);
            CHECK(db->standardTemperature() == 298.15);
            CHECK(db->elements().size() == 0);
            CHECK(db->species().size() == 0);
        }
    }

    if (flags & AIR5) {
        SECTION("5-species air") {
            db->load(std::string("N O N2 NO O2"));
            CHECK(db->standardPressure() == ONEATM);
            CHECK(db->standardTemperature() == 298.15);

            CHECK(db->elements().size() == 2);
            CHECK(db->elements()[0].name() == "N");
            CHECK(db->elements()[1].name() == "O");

            CHECK(db->species().size() == 5);
            CHECK(db->species()[0].name() == "N");
            CHECK(db->species()[1].name() == "O");
            CHECK(db->species()[2].name() == "N2");
            CHECK(db->species()[3].name() == "NO");
            CHECK(db->species()[4].name() == "O2");

            CHECK(db->species()[0].charge() == 0);
            CHECK(db->species()[0].molecularWeight() == Approx(mw_N));
            CHECK(db->species()[0].nAtoms() == 1);
            CHECK(db->species()[0].phase() == GAS);

            CHECK(db->species()[3].charge() == 0);
            CHECK(db->species()[3].molecularWeight() == Approx(mw_N + mw_O));
            CHECK(db->species()[3].nAtoms() == 2);
            CHECK(db->species()[3].phase() == GAS);
        }
    }

    if (flags & MULTIPHASE) {
        SECTION("Multiphase CO2") {
            db->load(std::string("C(gr) C O CO CO2"));
            CHECK(db->standardPressure() == ONEATM);
            CHECK(db->standardTemperature() == 298.15);

            CHECK(db->elements().size() == 2);
            CHECK(db->elements()[0].name() == "C");
            CHECK(db->elements()[1].name() == "O");

            CHECK(db->species().size() == 5);
            CHECK(db->species()[0].name() == "C");
            CHECK(db->species()[1].name() == "O");
            CHECK(db->species()[2].name() == "CO");
            CHECK(db->species()[3].name() == "CO2");
            CHECK(db->species()[4].name() == "C(gr)"); // solid phase last

            CHECK(db->species()[0].charge() == 0);
            CHECK(db->species()[0].molecularWeight() == Approx(mw_C));
            CHECK(db->species()[0].nAtoms() == 1);
            CHECK(db->species()[0].phase() == GAS);

            CHECK(db->species()[4].charge() == 0);
            CHECK(db->species()[4].molecularWeight() == Approx(mw_C));
            CHECK(db->species()[4].nAtoms() == 1);
            CHECK(db->species()[4].phase() == SOLID);
        }
    }
}

/**
 * Checks that data is loaded correctly from the NASA-7 thermodynamic database.
 */
TEST_CASE("Loading NASA-7 thermodynamic database", "[thermodynamics][loading]")
{
    ThermoDB* db = Utilities::Config::Factory<ThermoDB>::create("NASA-7", 0);
    checkThermoDBLoad(db, DEFAULTS | AIR5 | MULTIPHASE);
}

/**
 * Checks that data is loaded correctly from the NASA-9 thermodynamic database.
 */
TEST_CASE("Loading NASA-9 thermodynamic database", "[thermodynamics][loading]")
{
    ThermoDB* db = Utilities::Config::Factory<ThermoDB>::create("NASA-9", 0);
    checkThermoDBLoad(db, DEFAULTS | AIR5 | MULTIPHASE);
}

/**
 * Checks that data is loaded correctly from the new NASA-9 thermodynamic
 * database.
 */
TEST_CASE("Loading NASA-9-New thermodynamic database",
        "[thermodynamics][loading]"
)
{
    ThermoDB* db = Utilities::Config::Factory<ThermoDB>::create("NASA-9-New", 0);
    checkThermoDBLoad(db, DEFAULTS | AIR5);
}

/**
 * Checks that data is loaded correctly from the RRHO thermodynamic database.
 */
TEST_CASE("Loading RRHO thermodynamic database", "[thermodynamics][loading]")
{
    ThermoDB* db = Utilities::Config::Factory<ThermoDB>::create("RRHO", 0);
    checkThermoDBLoad(db, DEFAULTS | AIR5);
}


