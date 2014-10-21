/**
 * @file test_stefan_maxwell.cpp
 *
 * General tests on the Stefan Maxwell solver.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Stefan-Maxwell tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "MixtureFixtures.h"

// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(stefan_maxwell, MixtureFromCommandLine)

/*
 * Species diffusion fluxes must sum to zero.  This assumption depends on the
 * fact that the supplied driving forces also sum to zero so this must be done.
 */
BOOST_AUTO_TEST_CASE(rhoiV_sum_is_zero)
{
    std::vector<double> dp(mix().nSpecies());
    std::vector<double> Ji(mix().nSpecies());

    const double tol = 10.0*std::numeric_limits<double>::epsilon();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state
            equilibrateMixture(T, P);

            // Compute the driving forces for a temperature gradient in an
            // equilibrium mixture
            mix().dXidT(&dp[0]);

            // Compute the diffusion velocities
            double E; mix().stefanMaxwell(&dp[0], &Ji[0], E);
            double rho = mix().density();
            for (int i = 0; i < mix().nSpecies(); ++i)
                Ji[i] *= rho * mix().Y()[i];

            //for (int i = 0; i < mix().nSpecies(); ++i)
            //    std::cout << Ji[i] << std::endl;

            // Compute the max
            double max = std::abs(Ji[0]);
            for (int i = 1; i < mix().nSpecies(); ++i)
                max = std::max(std::abs(max), Ji[i]);

            // Make sure the sum is zero
            BOOST_CHECK_SMALL(
                std::accumulate(Ji.begin(), Ji.end(), 0.0)/max, tol);
        }
    }
}

// End of the stefan_maxwell test suite
BOOST_AUTO_TEST_SUITE_END()

