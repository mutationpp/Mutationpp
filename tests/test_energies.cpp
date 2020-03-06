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


TEST_CASE("Sum of species energies equals mixture energies",
    "[thermodynamics]"
)
{
    const double tol = std::numeric_limits<double>::epsilon()*1000;
    const int e_var_set = 0;
    const int t_var_set = 1;

    MIXTURE_LOOP
    (
        const int ns = mix.nSpecies();
        const int nt = mix.nEnergyEqns();

        std::vector<double> tmps1(nt);
        std::vector<double> tmps2(nt);

        Eigen::ArrayXd rhoi(ns);
        Eigen::ArrayXd species_energies(ns*nt);
        Eigen::ArrayXd mixture_energies(nt);

        EQUILIBRATE_LOOP
        (
            // Get the species densities
            double rho = mix.density();
            rhoi = rho * Eigen::Map<const Eigen::ArrayXd>(mix.Y(), ns);

            // Get a nominal temperature vector
            tmps1.assign(nt, T);

            // Compute a randomly perturbed temperature vector +- 500K
            for (int i = 0; i < nt; ++i)
                tmps1[i] += double(rand()) / RAND_MAX * 1000.0 - 500.0;

            // Set the state with the densities and perturbed temperatures
            mix.setState(rhoi.data(), &tmps1[0], t_var_set);

            // Check
            mix.getEnergiesMass(species_energies.data());
            mix.mixtureEnergies(mixture_energies.data());

            for (int i = 0; i < nt; ++i)
                CHECK(rho * mixture_energies[i] ==
                    Approx((rhoi*species_energies.segment(i*ns, ns)).sum()));
        )
    )
}


