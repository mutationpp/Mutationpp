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

/*
 * Checks that the thermal diffusion ratios sum to zero correctly such that
 * \f[
 * \sum_{i\in\mathcal{H}} k_{Ti} + k{Te}\frac{T_e}{T_h} = 0
 * \f]
 * in and out of thermal equilibrium.
 */
TEST_CASE("Thermal diffusion ratios sum to zero",
    "[transport]"
)
{
    const double tol = std::numeric_limits<double>::epsilon()*1000;

    MIXTURE_LOOP
    (
        VectorXd rhoi(mix.nSpecies());
        VectorXd tmps(mix.nEnergyEqns());
        VectorXd chii(mix.nSpecies());

        rhoi.setConstant(0.1);

        SECTION("Thermal equilibrium") {
            for (int i = 0; i < 20; ++i) {
                tmps.setConstant(500.0*i + 500.0);
                mix.setState(rhoi.data(), tmps.data(), 1);
                mix.heavyThermalDiffusionRatios(chii.data());
                CHECK(chii.sum() == Approx(0.0).margin(tol));
            }
        }

        SECTION("Thermal nonequilibrium") {
            for (int i = 0; i < 20; ++i) {
                double T = 500.0*i + 500.0;
                for (int k = 0; k < mix.nEnergyEqns(); ++k)
                    tmps[k] = T + 100.0*k*std::pow(-1.0,(double) k);
                mix.setState(rhoi.data(), tmps.data(), 1);
                mix.heavyThermalDiffusionRatios(chii.data());
                INFO("chii = " << chii);
                CHECK(chii.sum() == Approx(0.0).margin(tol));
            }
        }
    )
}


