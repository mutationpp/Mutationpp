/**
 * @file test_energies.cpp
 *
 * @brief General tests on the functions getting the energies.
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
#define BOOST_TEST_MODULE getEnergies() tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/random.hpp>

#include <Eigen/Dense>

#include "MixtureFixtures.h"

// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(energies, MixtureFromCommandLine)

/**
 * Tests that the sum of the species energies equals the mixture energies.
 */
BOOST_AUTO_TEST_CASE(test_species_energy_sum)
{
    const int ns = mix().nSpecies();
    const int nt = mix().nEnergyEqns();
    const int e_var_set = 0;
    const int t_var_set = 1;

    const double tol = 1.0e3*std::numeric_limits<double>::epsilon();

    std::vector<double> tmps1(nt);
    std::vector<double> tmps2(nt);

    Eigen::ArrayXd rhoi(ns);
    Eigen::ArrayXd species_energies(ns*nt);
    Eigen::ArrayXd mixture_energies(nt);

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
            double rho = mix().density();

            // Get the species densities
            rhoi = rho * Eigen::Map<const Eigen::ArrayXd>(mix().Y(), ns);

            // Get a nominal temperature vector
            tmps1.assign(nt, T);

            // Compute a randomly perturbed temperature vector
            for (int i = 0; i < nt; ++i)
                tmps1[i] += rand_temp();

            // Set the state with the densities and perturbed temperatures
            mix().setState(rhoi.data(), &tmps1[0], t_var_set);

            // Check
            mix().getEnergiesMass(species_energies.data());
            mix().mixtureEnergies(mixture_energies.data());

            for (int i = 0; i < nt; ++i)
                BOOST_CHECK_CLOSE(rho * mixture_energies[i],
                    (rhoi*species_energies.segment(i*ns, (i+1)*ns)).sum(), tol);
        }
    }
}

// End of the energies test suite
BOOST_AUTO_TEST_SUITE_END()
