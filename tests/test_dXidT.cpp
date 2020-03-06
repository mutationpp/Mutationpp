/**
 * @file test_dXidT.cpp
 *
 * @brief General tests on the dXidT function.
 */

/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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


TEST_CASE("Equilibrium mole fractions derivatives sum to zero",
    "[equilibrium][thermodynamics]"
)
{
    const double tol = std::numeric_limits<double>::epsilon()*1000;

    MIXTURE_LOOP
    (
        VectorXd dx(mix.nSpecies());

        EQUILIBRATE_LOOP
        (
            mix.dXidT(dx.data());
            INFO("dx = " << dx);
            REQUIRE(dx.sum() == Approx(0.0).margin(tol));
        )
    )
}

