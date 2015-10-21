/**
 * @file test_thermal_diff_ratios.cpp
 *
 * @brief Contains tests related to the species thermal diffusion ratios.
 * A mixture name must be given in the command line arguments for the tests to
 * work and all tests are performed on this mixture.
 */

/*
 * Copyright 2014-2015 von Karman Institute for Fluid Dynamics (VKI)
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
#define BOOST_TEST_MODULE thermal diffusion ratios tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <numeric>
#include "MixtureFixtures.h"


// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(thermal_diffusion_ratios, MixtureFromCommandLine)

/*
 * Checks that the thermal diffusion ratios sum to zero correctly such that
 * \f[
 * \sum_{i\in\mathcal{H}} k_{Ti} + k{Te}\frac{T_e}{T_h} = 0
 * \f]
 * in and out of thermal equilibrium.
 */
BOOST_AUTO_TEST_CASE(kTi_sum_to_zero)
{
    std::vector<double> rhoi(mix().nSpecies());
    std::vector<double> tmps(mix().nEnergyEqns());

    const double tol = 10.0*std::numeric_limits<double>::epsilon()*
            mix().nSpecies();

    // Check conditions with all temperatures the same ranging from 500 to
    // 10,000K
    for (int i = 0; i < 20; ++i) {
        rhoi.assign(mix().nSpecies(), 0.1);
        tmps.assign(mix().nEnergyEqns(), 500.0*i + 500.0);
        mix().setState(&rhoi[0], &tmps[0], 1);

        mix().thermalDiffusionRatios(&rhoi[0]);

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
            tmps[k] = 500.0*i + 500.0 + 100.0*k*std::pow(-1.0,(double) k);
        mix().setState(&rhoi[0], &tmps[0], 1);

        mix().thermalDiffusionRatios(&rhoi[0]);
        if (mix().hasElectrons())
            rhoi[0] *= mix().Te() / mix().T();

        double max = std::abs(rhoi[0]);
        for (int i = 1; i < mix().nSpecies(); ++i)
            max = std::max(std::abs(max), rhoi[i]);

        BOOST_CHECK_SMALL(
            std::accumulate(rhoi.begin(), rhoi.end(), 0.0)/max, tol);
    }
}

// End of the kinetics test suite
BOOST_AUTO_TEST_SUITE_END()


