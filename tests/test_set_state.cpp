/**
 * @file test_set_state.cpp
 *
 * @brief General tests on the setState() function.
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
#define BOOST_TEST_MODULE setState() tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/random.hpp>

#include "MixtureFixtures.h"

// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(set_state, MixtureFromCommandLine)

/*
 * All setState() functions (which take {rhoi, rho*Em} and {rhoi, Tm} as
 * variable sets should be able to compute rho*Em given a Tm and vice versa.
 * This test makes sure that you can go from one to the other and get the same
 * result.
 */
BOOST_AUTO_TEST_CASE(compare_T_from_E_from_T)
{
    const int ns = mix().nSpecies();
    const int nt = mix().nEnergyEqns();
    const int e_var_set = 0;
    const int t_var_set = 1;

    // @TODO investigate why this tolerance needs to be so big
    const double tol = 1.0e-8;//1.0e3*std::numeric_limits<double>::epsilon();

    std::vector<double> rhoi(ns);
    std::vector<double> tmps1(nt);
    std::vector<double> tmps2(nt);

    // Setup random number generator
    boost::mt19937 rng;
    boost::uniform_real<> temp_range(-500.0, 500.0);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
        rand_temp(rng, temp_range);

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state to make sure species temperatures and
            // densities are possible
            mix().equilibrate(T, P);

            // Get the species densities
            double rho = mix().density();
            for (int i = 0; i < ns; ++i)
                rhoi[i] = mix().Y()[i] * rho;

            // Get a nominal temperature vector
            tmps1.assign(nt, T);

            // Compute a randomly perturbed temperature vector
            for (int i = 0; i < nt; ++i)
                tmps1[i] += rand_temp();


            // Set the state with the densities and perturbed temperatures
            mix().setState(&rhoi[0], &tmps1[0], t_var_set);

            // Check that explicitly set properties still match
            mix().getTemperatures(&tmps2[0]);
            for (int i = 0; i < nt; ++i)
                BOOST_CHECK_CLOSE(tmps1[i], tmps2[i], tol);
            BOOST_CHECK_CLOSE(rho, mix().density(), tol);

            // Get the energy vector
            mix().mixtureEnergies(&tmps2[0]);
            for (int i = 0; i < nt; ++i)
                tmps2[i] *= rho;

            // Set the state of the mixture based on the energies
            mix().setState(&rhoi[0], &tmps2[0], e_var_set);

            // Check that implicitly set properties still match
            mix().getTemperatures(&tmps2[0]);
            for (int i = 0; i < nt; ++i)
                BOOST_CHECK_CLOSE(tmps1[i], tmps2[i], tol);
            BOOST_CHECK_CLOSE(rho, mix().density(), tol);
        }
    }
}

// End of the set_state test suite
BOOST_AUTO_TEST_SUITE_END()
