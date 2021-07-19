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

#include <catch.hpp>
#include "mutation++.h"
#include "Configuration.h"

#include <memory>

using namespace Catch;
using namespace Mutation;

// Make sure mixture options are set correctly
TEST_CASE("Mixture options", "[loading][mixtures]")
{
    GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Default options are set when using default constructor") {
        MixtureOptions opts;
        CHECK(opts.compositions().size() == 0);
        CHECK(opts.getDefaultComposition() == -1);
        CHECK(opts.getMechanism() == "none");
        CHECK(opts.getSource() == "");
        CHECK(opts.getSpeciesDescriptor() == "");
        CHECK(opts.getStateModel() == "ChemNonEq1T");
        CHECK(opts.getThermalConductivityAlgorithm() == "Chapmann-Enskog_LDLT");
        CHECK(opts.getThermodynamicDatabase() == "RRHO");
        CHECK(opts.getViscosityAlgorithm() == "Chapmann-Enskog_LDLT");
        CHECK(opts.hasDefaultComposition() == false);
    }

    SECTION("Default options are set when loading empty mixture") {
        MixtureOptions opts("default");
        CHECK(opts.compositions().size() == 0);
        CHECK(opts.getDefaultComposition() == -1);
        CHECK(opts.getMechanism() == "none");
        CHECK(opts.getSource() == Utilities::databaseFileName("default", "mixtures"));
        CHECK(opts.getSpeciesDescriptor() == "");
        CHECK(opts.getStateModel() == "ChemNonEq1T");
        CHECK(opts.getThermalConductivityAlgorithm() == "Chapmann-Enskog_LDLT");
        CHECK(opts.getThermodynamicDatabase() == "RRHO");
        CHECK(opts.getViscosityAlgorithm() == "Chapmann-Enskog_LDLT");
        CHECK(opts.hasDefaultComposition() == false);
    }

    SECTION("air11_NASA-7_ChemNonEq1T") {
        MixtureOptions opts("air11_NASA-7_ChemNonEq1T");

        CHECK(opts.compositions().size() == 1);
        CHECK(opts.compositions()[0].name() == "air");
        CHECK(opts.compositions()[0].size() == 2);
        CHECK(opts.compositions()[0][0].name == "N");
        CHECK(opts.compositions()[0][0].fraction == 0.79);
        CHECK(opts.compositions()[0][1].name == "O");
        CHECK(opts.compositions()[0][1].fraction == 0.21);
        CHECK(opts.getDefaultComposition() == 0);

        CHECK(opts.getMechanism() == "air11_mech");
        CHECK(opts.getSource() ==
            Utilities::databaseFileName("air11_NASA-7_ChemNonEq1T", "mixtures"));
        CHECK(opts.getSpeciesDescriptor() == "e- N N+ O O+ NO N2 N2+ O2 O2+ NO+");
        CHECK(opts.getStateModel() == "ChemNonEq1T");
        CHECK(opts.getThermalConductivityAlgorithm() == "Chapmann-Enskog_LDLT");
        CHECK(opts.getThermodynamicDatabase() == "NASA-7");
        CHECK(opts.getViscosityAlgorithm() == "Chapmann-Enskog_LDLT");
        CHECK(opts.hasDefaultComposition() == true);
    }

    SECTION("air5_RRHO_ChemNonEq1T_CG") {
        MixtureOptions opts("air5_RRHO_ChemNonEq1T_CG");

        CHECK(opts.compositions().size() == 1);
        CHECK(opts.compositions()[0].name() == "air");
        CHECK(opts.compositions()[0].size() == 2);
        CHECK(opts.compositions()[0][0].name == "N");
        CHECK(opts.compositions()[0][0].fraction == 0.79);
        CHECK(opts.compositions()[0][1].name == "O");
        CHECK(opts.compositions()[0][1].fraction == 0.21);
        CHECK(opts.getDefaultComposition() == 0);

        CHECK(opts.getMechanism() == "air5_mech");
        CHECK(opts.getSource() ==
            Utilities::databaseFileName("air5_RRHO_ChemNonEq1T_CG", "mixtures"));
        CHECK(opts.getSpeciesDescriptor() == "N O NO N2 O2");
        CHECK(opts.getStateModel() == "ChemNonEq1T");
        CHECK(opts.getThermalConductivityAlgorithm() == "Chapmann-Enskog_CG");
        CHECK(opts.getThermodynamicDatabase() == "RRHO");
        CHECK(opts.getViscosityAlgorithm() == "Chapmann-Enskog_CG");
        CHECK(opts.hasDefaultComposition() == true);
    }
}

// Check that basic data was successfully loaded for the test mixtures, this is
// just to make sure that if another test fails, that we should know if it was
// because the mixture didn't load.
void checkTestMixture(const std::string& name, int ns, int ne, int nr, bool he)
{
    SECTION(name) {
        Mixture mix(name);
        CHECK(mix.nSpecies() == ns);
        CHECK(mix.nElements() == ne);
        CHECK(mix.nReactions() == nr);
        CHECK(mix.nHeavy() == (he ? ns - 1 : ns));
        CHECK(mix.hasElectrons() == he);
    }
}

TEST_CASE("Loading test mixtures", "[loading][mixtures]")
{
    GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    checkTestMixture("air11_NASA-7_ChemNonEq1T", 11, 3, 9, true);
    checkTestMixture("air11_NASA-9_ChemNonEq1T", 11, 3, 9, true);
    checkTestMixture("air11_RRHO_ChemNonEq1T_CG", 11, 3, 9, true);
    checkTestMixture("air11_RRHO_ChemNonEq1T_Gupta-Yos", 11, 3, 9, true);
    checkTestMixture("air11_RRHO_ChemNonEq1T_Wilke", 11, 3, 9, true);
    checkTestMixture("air11_RRHO_ChemNonEq1T", 11, 3, 9, true);
    checkTestMixture("air11_RRHO_ChemNonEqTTv", 11, 3, 9, true);
    checkTestMixture("air5_NASA-7_ChemNonEq1T", 5, 2, 4, false);
    checkTestMixture("air5_NASA-9_ChemNonEq1T", 5, 2, 4, false);
    checkTestMixture("air5_RRHO_ChemNonEq1T_CG", 5, 2, 4, false);
    checkTestMixture("air5_RRHO_ChemNonEq1T_Gupta-Yos", 5, 2, 4, false);
    checkTestMixture("air5_RRHO_ChemNonEq1T_Wilke", 5, 2, 4, false);
    checkTestMixture("air5_RRHO_ChemNonEq1T", 5, 2, 4, false);
    checkTestMixture("air5_RRHO_ChemNonEqTTv", 5, 2, 4, false);
    checkTestMixture("argon_CR_ChemNonEq1T", 33, 2, 1984, true);
    checkTestMixture("argon_CR_ChemNonEqTTv", 33, 2, 1984, true);
}

class CoutRedirect
{
public:
    CoutRedirect() :
        m_buffer(),
        m_old_buffer(std::cout.rdbuf(m_buffer.rdbuf()))
    { }

    ~CoutRedirect() { std::cout.rdbuf(m_old_buffer); }

    std::string str() { return m_buffer.str(); }
private:
    std::stringstream m_buffer;
    std::streambuf* m_old_buffer;
};

void checkLoadMixture(
    const std::string& mix_name, int ns, int ne, int nr)
{
    std::unique_ptr<Mixture> mix;
    CHECK_NOTHROW(mix.reset(new Mixture(mix_name)));

    CHECK(mix->nSpecies() == ns);
    CHECK(mix->nElements() == ne);
    CHECK(mix->nReactions() == nr);

    {
        CoutRedirect cout;
        CHECK_NOTHROW(mix->collisionDB().Q11ij());
        INFO("standard out =\n" << cout.str());
        CHECK(cout.str().find("Warning") == std::string::npos);
    }
}

TEST_CASE("Loading air_5 mixture", "[loading][mixtures]")
{
    checkLoadMixture("air_5", 5, 2, 5);
}

TEST_CASE("Loading air_11 mixture", "[loading][mixtures]")
{
    checkLoadMixture("air_11", 11, 3, 22);
}

TEST_CASE("Loading Mars_19 mixture", "[loading][mixtures]")
{
    checkLoadMixture("Mars_19", 19, 4, 27);
}

TEST_CASE("Loading tacot-air_35 mixture", "[loading][mixtures]")
{
    checkLoadMixture("tacot-air_35", 35, 4, 0);
}
