/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include <catch.hpp>
#include <Eigen/Dense>

using namespace Mutation;
using namespace Catch;
using namespace Eigen;

TEST_CASE("stefanMaxwell yields diffusion fluxes which sum to zero",
    "[transport]"
)
{
    const double tol = std::numeric_limits<double>::epsilon()*1000;

    MIXTURE_LOOP
    (
        VectorXd dp(mix.nSpecies());
        VectorXd Ji(mix.nSpecies());

        EQUILIBRATE_LOOP
        (
            mix.dXidT(dp.data());

            // Scale to make the temperature gradient larger and make sure they
            // sum to zero
            dp *= 1.0 / dp.array().abs().maxCoeff();
            dp[0] -= dp.sum();

            // Compute the diffusion velocities
            double E; mix.stefanMaxwell(dp.data(), Ji.data(), E);
            double rho = mix.density();
            Ji.array() *= rho * Map<const ArrayXd>(mix.Y(), mix.nSpecies());

            // Make sure the sum is zero
            CHECK(Ji.sum() == Approx(0.0).margin(tol));
        )
    )
}

TEST_CASE("stefanMaxwell yields zero net conduction current",
    "[transport]"
)
{
    const double tol = std::numeric_limits<double>::epsilon()*1000;

    MIXTURE_LOOP
    (
        VectorXd dp(mix.nSpecies());
        VectorXd Ji(mix.nSpecies());

        EQUILIBRATE_LOOP
        (
            mix.dXidT(dp.data());

            // Scale to make the temperature gradient larger and make sure they
            // sum to zero
            dp *= 1.0 / dp.array().abs().maxCoeff();
            dp[0] -= dp.sum();

            // Compute the diffusion velocities
            double E; mix.stefanMaxwell(dp.data(), Ji.data(), E);
            for (int i = 0; i < mix.nSpecies(); ++i)
                Ji[i] *= mix.X()[i]*mix.speciesCharge(i);

            // Make sure the sum is zero
            CHECK(Ji.sum() == Approx(0.0).margin(tol));
        )
    )
}

