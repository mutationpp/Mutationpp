/**
 * @file test_wdot.cpp
 *
 * @brief Contains tests related to the species production rates in Mutation++.
 * A mixture name must be given in the command line arguments for the tests to
 * work and all tests are performed on this mixture and its associated
 * mechanism.
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
#define BOOST_TEST_MODULE wdot tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <numeric>
#include "MixtureFixtures.h"


// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(wdot, MixtureFromCommandLine)

/*
 * Ensures that the sum of the species production rates is zero for several
 * conditions.  This fact must be true if the reaction mechanism satisfies
 * elemental mass balance (which is checked when the mechanism is loaded).  It
 * has nothing to do with what the reaction rates are or what the species
 * densities are.  No matter what values are provided for these quantities, this
 * rule should hold, thus if it does not under any condition, there is a bug
 * present in the code.
 */
BOOST_AUTO_TEST_CASE(wdot_sum_equals_zero)
{
    std::vector<double> rhoi(mix().nSpecies());
    std::vector<double> tmps(mix().nEnergyEqns());

    const double tol = 10.0*std::numeric_limits<double>::epsilon();

    // Check conditions with all temperatures the same ranging from 500 to
    // 10,000K
    for (int i = 0; i < 20; ++i) {
        rhoi.assign(mix().nSpecies(), 0.1);
        tmps.assign(mix().nEnergyEqns(), 500.0*i + 500.0);
        mix().setState(&rhoi[0], &tmps[0], 1);
        mix().netProductionRates(&rhoi[0]);

        double max = std::abs(rhoi[0]);
        for (int i = 1; i < mix().nSpecies(); ++i)
            max = std::max(std::abs(max), rhoi[i]);

        BOOST_CHECK_SMALL(
            std::accumulate(rhoi.begin(), rhoi.end(), 0.0)/max, tol);
    }

    // Check conditions with different temperatures centered around temperatures
    // from 500 to 10,000K (thermal nonequilibrium)
    for (int i = 0; i < 20; ++i) {
        rhoi.assign(mix().nSpecies(), 0.1);
        for (int k = 0; k < mix().nEnergyEqns(); ++k)
            tmps[k] = 500.0*i + 500.0 + 50.0*k*std::pow(-1.0,(double) k);
        mix().setState(&rhoi[0], &tmps[0], 1);
        mix().netProductionRates(&rhoi[0]);

        double max = std::abs(rhoi[0]);
        for (int i = 1; i < mix().nSpecies(); ++i)
            max = std::max(std::abs(max), rhoi[i]);

        BOOST_CHECK_SMALL(
            std::accumulate(rhoi.begin(), rhoi.end(), 0.0)/max, tol);
    }
}

/*
 * Ensures that production rate is zero when mixture is in thermochemical
 * equilibrium for a range of temperatures and pressures.
 */
BOOST_AUTO_TEST_CASE(wdot_is_zero_in_equil)
{
    std::vector<double> rhoi(mix().nSpecies());
    std::vector<double> tmps(mix().nEnergyEqns());

    // Tolerance must account for errors in equilibrium solution which are
    // propagated forward into the production rates
    const double tol = 1.0e-10;

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state
            equilibrateMixture(T, P);

            // First compute the forward rate coefficient in order to scale the
            // results
            mix().forwardRateCoefficients(&rhoi[0]);
            double max = std::abs(rhoi[0]);
            for (int i = 1; i < mix().nSpecies(); ++i)
                max = std::max(std::abs(max), rhoi[i]);

            // Now compute the net production rates and make sure that the
            // scaled values are sufficiently small
            mix().netProductionRates(&rhoi[0]);
            for (int i = 0; i < mix().nSpecies(); ++i)
                BOOST_CHECK_SMALL(rhoi[i]/max, tol);
        }
    }
}

// End of the kinetics test suite
BOOST_AUTO_TEST_SUITE_END()


