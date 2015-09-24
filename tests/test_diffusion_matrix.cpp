/**
 * @file test_diffusion_matrix.cpp
 *
 * @brief General tests on the Diffusion matrix provider.
 */

/*
 * Copyright 2015 von Karman Institute for Fluid Dynamics (VKI)
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
using namespace Eigen;

#include "MixtureFixtures.h"

// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(diffusion_matrix, MixtureFromCommandLine)

/*
 * Species diffusion fluxes must sum to zero.  This assumption depends on the
 * fact that the supplied driving forces also sum to zero so this must be done.
 */
BOOST_AUTO_TEST_CASE(rhoiV_sum_is_zero)
{
    const int ns = mix().nSpecies();
    VectorXd dx(ns);
    VectorXd Ji(ns);

    const double tol = 1.0e-9;//10.0*std::numeric_limits<double>::epsilon();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state and compute the driving forces for a
            // temperature gradient in an equilibrium mixture
            mix().equilibrate(T, P);
            mix().dXidT(dx.data());

            // Scale to make the temperature gradient larger and make sure they
            // sum to zero
            //dx *= 1000.0 / dx.array().abs().maxCoeff();
            //dx[0] -= dx.sum();

            // Compute the diffusion velocities
            Ji = -mix().diffusionMatrix() * dx; // <- diffusion velocities
            Ji.array() *= mix().density() * Map<const ArrayXd>(mix().Y(),ns);

            // Make sure the sum is zero
            BOOST_CHECK_SMALL(Ji.sum()/(Ji.norm()*ns), tol);
        }
    }
}

// End of the stefan_maxwell test suite
BOOST_AUTO_TEST_SUITE_END()

