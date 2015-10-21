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

#include <Eigen/Dense>

#include "MixtureFixtures.h"

// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(stefan_maxwell, MixtureFromCommandLine)

/*
 * Species diffusion fluxes must sum to zero.  This assumption depends on the
 * fact that the supplied driving forces also sum to zero so this must be done.
 */
BOOST_AUTO_TEST_CASE(rhoiV_sum_is_zero)
{
    Eigen::VectorXd dp(mix().nSpecies());
    Eigen::VectorXd Ji(mix().nSpecies());

    const double tol = 1.0e-9;//10.0*std::numeric_limits<double>::epsilon();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state and compute the driving forces for a
            // temperature graient in an equilibrium mixture
            mix().equilibrate(T, P);
            mix().dXidT(dp.data());

            // Scale to make the temperature gradient larger and make sure they
            // sum to zero
            dp *= 1000.0 / dp.array().abs().maxCoeff();
            dp[0] -= dp.sum();

            // Compute the diffusion velocities
            double E; mix().stefanMaxwell(&dp[0], &Ji[0], E);
            double rho = mix().density();
            for (int i = 0; i < mix().nSpecies(); ++i)
                Ji[i] *= rho * mix().Y()[i];

            // Compute the norm
            double norm = Ji.norm() * mix().nSpecies();

            // Make sure the sum is zero
            BOOST_CHECK_SMALL(Ji.sum()/norm, tol);
        }
    }
}

/*
 * There should be no net conduction current. Tests that sum of Xi*qi*Vi is
 * small.
 */
BOOST_AUTO_TEST_CASE(XiqiVi_sum_is_zero)
{
    Eigen::VectorXd dp(mix().nSpecies());
    Eigen::VectorXd Ji(mix().nSpecies());

    const double tol = 1.0e-9;//10.0*std::numeric_limits<double>::epsilon();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 15000.0;

            // Set an equilibrium state and compute the driving forces for a
            // temperature graient in an equilibrium mixture
            mix().equilibrate(T, P);
            mix().dXidT(dp.data());

            // Scale to make the temperature gradient larger and make sure they
            // sum to zero
            dp *= 1000.0 / dp.array().abs().maxCoeff();
            dp[0] -= dp.sum();

            // Compute the charge fluxes and norm of the velocities
            double E; mix().stefanMaxwell(&dp[0], &Ji[0], E);
            double norm = Ji.norm()*Mutation::QE;
            for (int i = 0; i < mix().nSpecies(); ++i)
                Ji[i] *= mix().X()[i]*mix().speciesCharge(i);

            // Make sure the sum is zero
            BOOST_CHECK_SMALL(Ji.sum()/norm, tol);
        }
    }
}

// End of the stefan_maxwell test suite
BOOST_AUTO_TEST_SUITE_END()

