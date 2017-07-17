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

#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include "catch.hpp"
#include <Eigen/Dense>

using namespace Mutation;
using namespace Catch;
using namespace Eigen;

/*
 * Ensures that the sum of the species production rates is zero for several
 * conditions.  This fact must be true if the reaction mechanism satisfies
 * elemental mass balance (which is checked when the mechanism is loaded).  It
 * has nothing to do with what the reaction rates are or what the species
 * densities are.  No matter what values are provided for these quantities, this
 * rule should hold, thus if it does not under any condition, there is a bug
 * present in the code.
 */
TEST_CASE
(
    "Species production rates sum to zero",
    "[kinetics]"
)
{
    const double tol = std::numeric_limits<double>::epsilon();

    MIXTURE_LOOP
    (
        VectorXd rhoi(mix.nSpecies());
        VectorXd tmps(mix.nEnergyEqns());

        SECTION("Thermal equilibrium") {
            for (int i = 0; i < 20; ++i) {
                rhoi.setConstant(0.1);
                tmps.setConstant(500.0*i + 500.0);
                mix.setState(rhoi.data(), tmps.data(), 1);
                mix.netProductionRates(rhoi.data());
                INFO("rhoi = " << rhoi);
                REQUIRE(rhoi.sum() == Approx(0.0).margin(tol));
            }
        }

        SECTION("Thermal nonequilibrium") {
            for (int i = 0; i < 20; ++i) {
                rhoi.setConstant(0.1);
                for (int k = 0; k < mix.nEnergyEqns(); ++k)
                    tmps[k] = 500.0*i + 500.0 + 50.0*k*std::pow(-1.0,(double) k);
                mix.setState(rhoi.data(), tmps.data(), 1);
                mix.netProductionRates(rhoi.data());
                CHECK(rhoi.sum() == Approx(0.0).margin(tol));
            }
        }
    )
}



/*
 * Ensures that production rate is zero when mixture is in thermochemical
 * equilibrium for a range of temperatures and pressures.
 */
//BOOST_AUTO_TEST_CASE(wdot_is_zero_in_equil)
//{
//    std::vector<double> rhoi(mix().nSpecies());
//    std::vector<double> tmps(mix().nEnergyEqns());
//
//    // Tolerance must account for errors in equilibrium solution which are
//    // propagated forward into the production rates
//    const double tol = 1.0e-7;
//
//    // Loop over range of pressures and temperatures
//    for (int ip = 0.0; ip < 10; ++ip) {
//        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
//        for (int it = 0.0; it < 10; ++it) {
//            double T = 1000.0*it + 1000.0;
//
//            // Set an equilibrium state
//            mix().equilibrate(T, P);
//
//            mix().densities(&rhoi[0]);
//            double max = rhoi[0];
//            for (int i = 1; i < mix().nSpecies(); ++i)
//                max = std::max(max, rhoi[i]);
//
//            // Now compute the net production rates and make sure that the
//            // scaled values are sufficiently small
//            mix().netProductionRates(&rhoi[0]);
//            for (int i = 0; i < mix().nSpecies(); ++i)
//                BOOST_CHECK_SMALL(rhoi[i]/max, tol);
//        }
//    }
//}
//
//// End of the kinetics test suite
//BOOST_AUTO_TEST_SUITE_END()


