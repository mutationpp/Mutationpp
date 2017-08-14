/*
 * Copyright 2017 von Karman Institute for Fluid Dynamics (VKI)
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
#include <catch/catch.hpp>
#include <eigen3/Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Utilities;
using namespace Catch;
using namespace Eigen;

TEST_CASE
(
    "Units",
    "[utilities]"
)
{
    Units units;
    CHECK_NOTHROW(units = Units("kg -m/s-s"));
    CHECK_NOTHROW(units = Units(std::string("kg-m /s-s")));
    CHECK_NOTHROW(units = Units(Units("kg-m/s - s")));
    CHECK_THROWS_AS(Units("text"), InvalidInputError);

    try {
        units = Units("kg-m/b-b");
    } catch (InvalidInputError& e) {
        CHECK(e.inputName() == "units");
        CHECK(e.inputValue() == "kg-m/b-b");
    }

    CHECK_NOTHROW(units = Units("g/cm-cm-cm"));
    CHECK(units.convertToBase(1.0) == Approx(1000000.0 / 1000.0));
    CHECK(units.convertTo(1.0, "kg/m-m-m") == units.convertToBase(1.0));
    CHECK(units.convertTo(1.0, units) == 1.0);

    CHECK(Units("kg-m/s-s").simplestRepresentation() == "N");
    CHECK(Units("N-m").simplestRepresentation() == "J");

    CHECK(units.convertTo(1.0, "cal") == 1.0);
}
