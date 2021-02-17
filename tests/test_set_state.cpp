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
 * All setState() functions (which take {rhoi, rho*Em} and {rhoi, Tm} as
 * variable sets should be able to compute rho*Em given a Tm and vice versa.
 * This test makes sure that you can go from one to the other and get the same
 * result.
 */
TEST_CASE("setState converts rho*Em to Tm and vice a versa",
        "[thermodynamics]"
)
{
    const int e_var_set = 0;
    const int t_var_set = 1;

    MIXTURE_LOOP
    (
        const int ns = mix.nSpecies();
        const int nt = mix.nEnergyEqns();

        VectorXd rhoi(ns);
        VectorXd T1(nt);
        VectorXd T2(nt);
        VectorXd e1(nt);
        VectorXd e2(nt);

        EQUILIBRATE_LOOP
        (
            // Get the species densities
            double rho = mix.density();
            rhoi = rho * Eigen::Map<const Eigen::ArrayXd>(mix.Y(), ns);

            // Compute a randomly perturbed temperature vector T +- 500K
            for (int i = 0; i < nt; ++i)
                T1[i] = double(rand()) / RAND_MAX * 1000.0 - 500.0 + T;

            // Set the state with the densities and perturbed temperatures
            mix.setState(rhoi.data(), T1.data(), t_var_set);

            // Check that explicitly set properties still match
            mix.getTemperatures(T2.data());
            for (int i = 0; i < nt; ++i)
                CHECK(T1[i] == Approx(T2[i]));
            
            CHECK(rho == Approx(mix.density()));

            // Get the energy vector
            mix.mixtureEnergies(e1.data());
            e1 *= rho;

            // Set the state of the mixture based on the energies
            mix.setState(rhoi.data(), e1.data(), e_var_set);

            // Check that implicitly set properties still match
            mix.getTemperatures(T2.data());
            for (int i = 0; i < nt; ++i)
                CHECK(T2[i] ==  Approx(T2[i]));
            
            CHECK(rho == Approx(mix.density()));

            mix.mixtureEnergies(e2.data());
            e2 *= rho;
            for (int i = 0; i < nt; ++i)
                CHECK(e1[i] ==  Approx(e2[i]));
        )
    )
}

