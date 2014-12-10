/**
 * @file test_stefan_maxwell.cpp
 *
 * @brief General tests on the Stefan Maxwell solver.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

    const double tol = 1.0e-8;//10.0*std::numeric_limits<double>::epsilon();

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
            // Scale to make the temperature gradient larger
            double max = std::abs(dp[0]);
            for (int i = 1; i < mix().nSpecies(); ++i)
                max = std::max(max, std::abs(dp[i]));
            for (int i = 0; i < mix().nSpecies(); ++i)
                dp[i] *= 1000.0/max;

            // Compute the diffusion velocities
            double E; mix().stefanMaxwell(&dp[0], &Ji[0], E);
            double rho = mix().density();
            for (int i = 0; i < mix().nSpecies(); ++i)
                Ji[i] *= rho * mix().Y()[i];

            //for (int i = 0; i < mix().nSpecies(); ++i)
            //    std::cout << Ji[i] << std::endl;

            // Compute the norm
            double norm = 0.0;
            for (int i = 0; i < mix().nSpecies(); ++i)
                norm += Ji[i]*Ji[i];
            norm = std::sqrt(norm);

            // Make sure the sum is zero
            BOOST_CHECK_SMALL(
                    std::accumulate(Ji.begin(), Ji.end(), 0.0)/norm, tol);
        }
    }
}

/*
 * Electron diffusion velocity must be zero (only tested if electrons are
 * present in the mixture.
 */
//BOOST_AUTO_TEST_CASE(Ve_is_zero)
//{
//    if (!mix().hasElectrons())
//        return;
//
//    std::vector<double> dp(mix().nSpecies());
//    std::vector<double> Vi(mix().nSpecies());
//
//    const double tol = 10.0*std::numeric_limits<double>::epsilon();
//
//    // Loop over range of pressures and temperatures
//    for (int ip = 0.0; ip < 10; ++ip) {
//        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
//        for (int it = 0.0; it < 10; ++it) {
//            double T = 1000.0*it + 1000.0;
//
//            // Set an equilibrium state
//            equilibrateMixture(T, P);
//
//            // Compute the driving forces for a temperature gradient in an
//            // equilibrium mixture
//            mix().dXidT(&dp[0]);
//
//            // Compute the diffusion velocities
//            double E; mix().stefanMaxwell(&dp[0], &Vi[0], E);
//
//            // Make sure Ve is zero
//            BOOST_CHECK_SMALL(Vi[0], tol);
//        }
//    }
//}

/*
 * There should be no net conduction current. Tests that sum of Xi*qi*Vi is
 * small.
 */
BOOST_AUTO_TEST_CASE(XiqiVi_sum_is_zero)
{
    std::vector<double> dp(mix().nSpecies());
    std::vector<double> Ji(mix().nSpecies());

    const double tol = 1.0e-9;//10.0*std::numeric_limits<double>::epsilon();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 15000.0;

            // Set an equilibrium state
            equilibrateMixture(T, P);

            // Compute the driving forces for a temperature gradient in an
            // equilibrium mixture
            mix().dXidT(&dp[0]);
            // Scale to make the temperature gradient larger
            double max = std::abs(dp[0]);
            for (int i = 1; i < mix().nSpecies(); ++i)
                max = std::max(max, std::abs(dp[i]));
            for (int i = 0; i < mix().nSpecies(); ++i)
                dp[i] *= 1000.0*max;

            // Compute the charge fluxes and norm of the velocities
            double E; mix().stefanMaxwell(&dp[0], &Ji[0], E);
            double norm = 0.0;
            for (int i = 0; i < mix().nSpecies(); ++i) {
                norm += Ji[i]*Ji[i];
                Ji[i] *= mix().X()[i]*mix().speciesCharge(i);
            }
            norm = std::sqrt(norm)*Mutation::QE;

            // Make sure the sum is zero
            BOOST_CHECK_SMALL(
                    std::accumulate(Ji.begin(), Ji.end(), 0.0)/norm, tol);
        }
    }
}

// End of the stefan_maxwell test suite
BOOST_AUTO_TEST_SUITE_END()

