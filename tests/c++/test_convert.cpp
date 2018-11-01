/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
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


TEST_CASE("Converting thermodynamic properties of species and elements",
        "[thermodynamic]"
)
{
    const double tol = std::numeric_limits<double>::epsilon();

    MIXTURE_LOOP
    (
        Eigen::VectorXd Ys(mix.nSpecies());
        Eigen::VectorXd rhoi(mix.nSpecies());

        Eigen::VectorXd Ye(mix.nElements());
        Eigen::VectorXd Xe(mix.nElements());

        EQUILIBRATE_LOOP
        (
            mix.densities(rhoi.data());

            mix.convert<Mutation::Thermodynamics::RHO_TO_CONC>(rhoi.data(),Ys.data());

            mix.convert<Mutation::Thermodynamics::CONC_TO_Y>(Ys.data(),Ys.data());

            CHECK(Ys.sum() == Approx(1.0).margin(tol));

            mix.convert<Mutation::Thermodynamics::Y_TO_X>(Ys.data(),Ys.data());

            CHECK(Ys.sum() == Approx(1.0).margin(tol));

            mix.convert<Mutation::Thermodynamics::X_TO_Y>(Ys.data(),Ys.data());

            CHECK(Ys.sum() == Approx(1.0).margin(tol));

            mix.convert<Mutation::Thermodynamics::Y_TO_YE>(Ys.data(),Ye.data());

            CHECK(Ys.sum() == Approx(1.0).margin(tol));

            mix.convert<Mutation::Thermodynamics::YE_TO_XE>(Ye.data(),Xe.data());

            CHECK(Ys.sum() == Approx(1.0).margin(tol));

            CHECK( Xe.isApprox(Eigen::Map<const Eigen::VectorXd>(mix.getDefaultComposition(),mix.nElements())));
        )
    )
}
