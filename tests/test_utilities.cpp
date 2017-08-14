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
#include <catch/catch.hpp>
#include <eigen3/Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Utilities;
using namespace Catch;
using namespace Eigen;

typedef Matrix<double,7,1> Dimension;

void checkUnits(const Units& units, Dimension dimension, double factor = 1.0)
{
    INFO("factor = " << factor);
    INFO("dimension = " << dimension.transpose());
    INFO("exponents = " << units.exponents().transpose());
    CHECK(units.factor() == Approx(factor));
    CHECK(units.exponents().isApprox(dimension));
}

void checkUnits(
    const char* const symbol, Dimension dimension, double factor = 1.0)
{
    INFO("symbol = " << symbol);
    INFO("factor = " << factor);
    INFO("dimension = " << dimension.transpose());
    Units units;
    CHECK_NOTHROW(units = Units(symbol));
    checkUnits(units, dimension, factor);
}

/**
 * Tests the basic functionality of the Units type.
 */
TEST_CASE
(
    "Units",
    "[utilities]"
)
{
    // Define the expected dimensions for different types of units
    Dimension length      = MatrixXd::Identity(7,7).col(0);
    Dimension mass        = MatrixXd::Identity(7,7).col(1);
    Dimension time        = MatrixXd::Identity(7,7).col(2);
    Dimension current     = MatrixXd::Identity(7,7).col(3);
    Dimension temperature = MatrixXd::Identity(7,7).col(4);
    Dimension quantity    = MatrixXd::Identity(7,7).col(5);
    Dimension intensity   = MatrixXd::Identity(7,7).col(6);
    Dimension force       = length + mass - 2*time;
    Dimension energy      = force + length;
    Dimension power       = energy - time;
    Dimension area        = 2*length;
    Dimension volume      = 3*length;
    Dimension density     = mass - volume;
    Dimension pressure    = force - area;

    // Check them dimensions and factors of predefined units
    checkUnits("Å"       , length, 1.0e-10);
    checkUnits("mm"      , length, 0.001);
    checkUnits("cm"      , length, 0.01);
    checkUnits("m"       , length);
    checkUnits("km"      , length, 1000.0);
    checkUnits("g"       , mass, 0.001);
    checkUnits("kg"      , mass);
    checkUnits("s"       , time);
    checkUnits("min"     , time, 60.0);
    checkUnits("A"       , current);
    checkUnits("K"       , temperature);
    checkUnits("molecule", quantity, 6.0221415e-23);
    checkUnits("mol"     , quantity);
    checkUnits("kmol"    , quantity, 1000.0);
    checkUnits("cd"      , intensity);
    checkUnits("N"       , force);
    checkUnits("meV"     , energy, 6.242e-21);
    checkUnits("eV"      , energy, 6.242e-18);
    checkUnits("J"       , energy);
    checkUnits("cal"     , energy, 4.184);
    checkUnits("kJ"      , energy, 1000.0);
    checkUnits("kcal"    , energy, 4184.0);
    checkUnits("Pa"      , pressure);
    checkUnits("bar"     , pressure, 100000.0);
    checkUnits("atm"     , pressure, 101325.0);

    // Check some combinations of units
    checkUnits("kg m/s s"  , force);
    checkUnits("cal / s"   , power, 4.184);
    checkUnits("g/cm-cm-cm", density, 1000.0);
    checkUnits("km-Å"      , area, 1.0e-7);
    checkUnits("cm.K"      , length + temperature, 0.01);
    checkUnits(Units("cm")^3, volume, 1.0e-6);
    checkUnits(Units("kg")*Units("m")/(Units("s")^2)*100.0, force, 100.0);

    // Check that invalid input fails correctly
    CHECK_THROWS_AS(Units(" "   )   , InvalidInputError);
    CHECK_THROWS_AS(Units("/"   )   , InvalidInputError);
    CHECK_THROWS_AS(Units(" /m-m")  , InvalidInputError);
    CHECK_THROWS_AS(Units("m-m/ ")  , InvalidInputError);
    CHECK_THROWS_AS(Units("m/kg/s") , InvalidInputError);
    CHECK_THROWS_AS(Units("text")   , InvalidInputError);
    CHECK_THROWS_AS(Units("m^2" )   , InvalidInputError);

    // Check conversions
    CHECK(Units("cal").convertTo(1.0, "kJ") == Approx(4.184 / 1000.0));
    CHECK(Units("Å.cm/s").convertTo(1.0, "m.m/s") == Approx(1.0e-12));
    CHECK(Units("g/cm.cm.cm").convertToBase(1.0) == Approx(1000.0));

    // String conversion
    CHECK(Units("kg-m/s-s").toString() == "m kg / s^2");
    CHECK(Units("N-m").toString() == "m^2 kg / s^2");
}
