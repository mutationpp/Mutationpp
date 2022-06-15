#include "MillikanWhite.h"
#include "Mixture.h"
#include "TemporaryFile.h"
#include <catch.hpp>

using namespace Mutation;
using namespace Mutation::Transfer;
using namespace Catch;

void checkDefaultModelData(bool with_electrons)
{
    const std::string species_list = (with_electrons ? "e- N2 O2 N" : "N2 O2 N");
    const size_t offset = (with_electrons ? 1 : 0);
        
    MixtureOptions opts;
    opts.setSpeciesDescriptor(species_list);
    Mixture mix(opts);

    const double thetav = 3408.464;
    MillikanWhiteModelData data(mix, offset, thetav);
    
    CHECK(data.speciesIndex() == offset);
    CHECK(data.molecularWeight() == mix.speciesMw(offset));
    CHECK(data.referenceCrossSection() == 1.0e-20);
    CHECK(data.limitingCrossSection(mix.T()) == 1.0e-20 * (2.5e9/(mix.T()*mix.T())));
    REQUIRE(data.a().size() == 3);
    REQUIRE(data.b().size() == 3);

    const double mw = mix.speciesMw(offset);

    for (size_t i = 0; i < 3; ++i)
    {
        const double mwi = mix.speciesMw(i+offset);
        const double mu = 1000.0*mw*mwi/(mw+mwi); // g/mol
        
        CHECK(data.a()[i] == 1.16e-3*std::sqrt(mu)*std::pow(thetav, 4.0/3.0));
        CHECK(data.b()[i] == 0.015*std::pow(mu, 0.25));
    }
}

/**
 * Tests that Millikan-White model data defaults to correct values.
 */
TEST_CASE("MillikanWhiteModelData provides default model values", "[transfer]")
{
    checkDefaultModelData(false);
    checkDefaultModelData(true);
}

void checkDefaultRelaxationRate(bool with_electrons)
{
    INFO(std::string("with_electrons = ") + (with_electrons ? "true" : "false"));
    const std::string species_list = (with_electrons ? "e- N2" : "N2");
    const size_t offset = (with_electrons ? 1 : 0);
    
    MixtureOptions opts;
    opts.setSpeciesDescriptor(species_list);
    opts.setStateModel("ChemNonEq1T");
    Mixture mix(opts);

    const double thetav = 3408.464;
    MillikanWhiteModel model({mix, offset, thetav});

    const int SET_Y_AND_PT = 2;
    
    std::vector<double> y;
    if (with_electrons)
        y.push_back(0.0);
    y.push_back(1.0);

    const double pT[] = { ONEATM, 300.0 };
    mix.setState(y.data(), pT, SET_Y_AND_PT);
    
    const double mu = 1000.0*0.5*mix.speciesMw(offset);
    const double A = 1.16e-3*std::sqrt(mu)*std::pow(thetav, 4.0/3.0);
    const double B = 0.015*std::pow(mu, 0.25);
    const double tau_mw = std::exp(A*(std::pow(300.0, -1.0/3.0) - B) - 18.42);
    
    const double ni = mix.numberDensity();
    const double ci = std::sqrt(8*RU*300.0/(PI*mix.speciesMw(offset)));
    const double sigma = 1e-20 * (2.5E9/(300.0*300.0));
    const double tau_park = 1.0/(ni * ci * sigma);
    
    CHECK(model.relaxationTime(mix) == tau_mw + tau_park);
}

TEST_CASE("MillikanWhiteModel provides correct default relaxation time", "[transfer]")
{
    checkDefaultRelaxationRate(false);
    checkDefaultRelaxationRate(true);
}

TEST_CASE("Can load MillikanWhiteModel from database", "[transfer]")
{
    // Create a temporary database
    Mutation::Utilities::IO::TemporaryFile file(".xml");
    file << R"(
    <VT>
        <Millikan-White>
            <vibrator species="N2" omegav="3.0E-21">
                <partner species="N" a="180.0" b ="0.0262" />
                <partner species="N2" a="221" b="0.0290" />
            </vibrator>
        </Millikan-White>
    </VT>)";
    file.close();

    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2 N");
    Mixture mix(opts);

    MillikanWhiteModelDB db(mix, file.filename());

    auto model = db.create("N2", 3408.464);

    CHECK(model.speciesIndex() == 0);
    CHECK(model.molecularWeight() == mix.speciesMw(0));
    CHECK(model.data().a().size() == 2);
    CHECK(model.data().b().size() == 2);

    CHECK(model.data().a()[0] == 221.0);
    CHECK(model.data().b()[0] == 0.0290);
    
    CHECK(model.data().a()[1] == 180.0);
    CHECK(model.data().b()[1] == 0.0262);
}

