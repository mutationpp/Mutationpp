#include "HarmonicOscillator.h"
#include "TemporaryFile.h"
#include <catch.hpp>

using namespace Mutation::Thermodynamics;
using namespace Catch;

/**
 * Tests construction of harmonic oscillators.
 */
TEST_CASE("Harmonic oscillator construction", "[thermo]")
{
    SECTION("No characteristic temperatures")
    {
        HarmonicOscillator ho;
        CHECK(ho.characteristicTemperatures().size() == 0);
    }
    
    SECTION("One characteristic temperature")
    {
        std::vector<double> temps = { 1.0 };
        HarmonicOscillator ho(temps);
        CHECK(ho.characteristicTemperatures().size() == 1);
        CHECK(ho.characteristicTemperatures()[0] == 1.0);
    }

    SECTION("One characteristic temperature")
    {
        std::vector<double> temps = { 1.0, 2.0, 3.0 };
        HarmonicOscillator ho(temps);
        CHECK(ho.characteristicTemperatures().size() == 3);
        CHECK(ho.characteristicTemperatures()[0] == 1.0);
        CHECK(ho.characteristicTemperatures()[1] == 2.0);
        CHECK(ho.characteristicTemperatures()[2] == 3.0);
    }
}

/**
 * Tests harmonic oscillators energy calculation.
 */
TEST_CASE("Harmonic oscillator energy", "[thermo]")
{
    SECTION("No characteristic temperatures")
    {
        HarmonicOscillator ho;
        CHECK(ho.energy(1.0) == 0.0);
    }
    
    SECTION("One characteristic temperature")
    {
        std::vector<double> temps = { 1.0 };
        HarmonicOscillator ho(temps);
        CHECK(ho.energy(1.0) == Approx(0.581976706869326));
    }

    SECTION("One characteristic temperature")
    {
        std::vector<double> temps = { 1.0, 2.0, 3.0 };
        HarmonicOscillator ho(temps);
        CHECK(ho.energy(1.0) == Approx(1.052199081842426));
    }
}


/**
 * Tests loading of harmonic oscillators from database.
 */
TEST_CASE("Loading harmonic oscillator from database", "[thermo]")
{
    // Create a temporary database
    Mutation::Utilities::IO::TemporaryFile file(".xml");
    file << R"(
    <specieslist>
        <species name="CO2">
            <thermodynamics type="RRHO">
                <vibrational_temperatures units="K">
                    932.109 932.109 1914.081 3373.804
                </vibrational_temperatures>
            </thermodynamics>
        </species>
        <species name="N2">
            <thermodynamics type="RRHO">
                <vibrational_temperatures units="K">
                    3408.464
                </vibrational_temperatures>
            </thermodynamics>
        </species>
    </specieslist>)";
    file.close();

    HarmonicOscillatorDB db(file.filename());

    auto ho = db.create("CO2");
    CHECK(ho.characteristicTemperatures().size() == 4);
    CHECK(ho.characteristicTemperatures()[0] == 932.109);
    CHECK(ho.characteristicTemperatures()[1] == 932.109);
    CHECK(ho.characteristicTemperatures()[2] == 1914.081);
    CHECK(ho.characteristicTemperatures()[3] == 3373.804);

    ho = db.create("N2");
    CHECK(ho.characteristicTemperatures().size() == 1);
    CHECK(ho.characteristicTemperatures()[0] == 3408.464);

    ho = db.create("test");
    CHECK(ho.characteristicTemperatures().size() == 0);
}


