/**
 * @file MixtureFixtures.h
 *
 * Provides test fixtures which provide preloaded mixture for different test
 * cases to use.
 */

#include <boost/test/included/unit_test.hpp>
#include "mutation++.h"

class MixtureFromCommandLine {
public:
    MixtureFromCommandLine()
        : mp_mix(NULL)
    {
        BOOST_REQUIRE_MESSAGE(
            boost::unit_test::framework::master_test_suite().argc > 1,
            "Mixture name must be given in command line arguments");

        // Load the mixture from the first command line argument
        BOOST_TEST_MESSAGE("Loading mixture from command line argument");
        mp_mix = new Mutation::Mixture(
            boost::unit_test::framework::master_test_suite().argv[1]);
    }

    ~MixtureFromCommandLine()
    {
        // Load the mixture from the first command line argument
        BOOST_TEST_MESSAGE("Destroying mixture from command line argument");
        delete mp_mix;
    }

    Mutation::Mixture& mix() { return *mp_mix; }

    /**
     * Sets the mixture to an equilibrium state at the given T and P.  Assumes
     * that the variable set 1 of the state model accepts species densities and
     * temperature vectors.
     */
    void equilibrateMixture(double T, double P)
    {
        static std::vector<double> rhoi(mix().nSpecies());
        static std::vector<double> tmps(mix().nEnergyEqns());

        // Get equilibrium mole fractions
        mix().equilibriumComposition(T, P, &rhoi[0]);

        // Convert to species densities
        double rho = mix().density(T, P, &rhoi[0]);
        mix().convert<Mutation::Thermodynamics::X_TO_Y>(&rhoi[0], &rhoi[0]);
        for (int i = 0; i < mix().nSpecies(); ++i)
            rhoi[i] *= rho;

        tmps.assign(mix().nEnergyEqns(), T);
        mix().setState(&rhoi[0], &tmps[0], 1);
    }

private:
    Mutation::Mixture* mp_mix;
};
