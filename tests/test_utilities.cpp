/*
 * Copyright 2017-2020 von Karman Institute for Fluid Dynamics (VKI)
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
#include <catch.hpp>
#include <Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
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
TEST_CASE("Units", "[utilities]")
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

/**
 * Tests for utility functions.
 */
TEST_CASE("Utility functions", "[utilities]")
{
    SECTION("ordinalSuffix()") {
        CHECK(ordinalSuffix(  0) == "th");
        CHECK(ordinalSuffix(  1) == "st");
        CHECK(ordinalSuffix(  2) == "nd");
        CHECK(ordinalSuffix(  3) == "rd");
        CHECK(ordinalSuffix(  4) == "th");
        CHECK(ordinalSuffix(  5) == "th");
        CHECK(ordinalSuffix(  6) == "th");
        CHECK(ordinalSuffix(  7) == "th");
        CHECK(ordinalSuffix(  8) == "th");
        CHECK(ordinalSuffix(  9) == "th");
        CHECK(ordinalSuffix( 10) == "th");
        CHECK(ordinalSuffix( 11) == "th");
        CHECK(ordinalSuffix( 12) == "th");
        CHECK(ordinalSuffix( 13) == "th");
        CHECK(ordinalSuffix( 14) == "th");
        CHECK(ordinalSuffix( 15) == "th");
        CHECK(ordinalSuffix( 16) == "th");
        CHECK(ordinalSuffix( 17) == "th");
        CHECK(ordinalSuffix( 18) == "th");
        CHECK(ordinalSuffix( 19) == "th");
        CHECK(ordinalSuffix( 20) == "th");
        CHECK(ordinalSuffix( 21) == "st");
        CHECK(ordinalSuffix( 22) == "nd");
        CHECK(ordinalSuffix( 23) == "rd");
        CHECK(ordinalSuffix(501) == "st");
        CHECK(ordinalSuffix(502) == "nd");
        CHECK(ordinalSuffix(503) == "rd");
        CHECK(ordinalSuffix(504) == "th");
        CHECK(ordinalSuffix(511) == "th");
        CHECK(ordinalSuffix(512) == "th");
        CHECK(ordinalSuffix(513) == "th");
        CHECK(ordinalSuffix(521) == "st");
        CHECK(ordinalSuffix(522) == "nd");
        CHECK(ordinalSuffix(523) == "rd");
    }
}


/**
 * Tests the XML classes
 */
TEST_CASE("XML classes", "[utilities]")
{
    TemporaryFile file;
    file << "<tag att=\"100\">\n"
         << "    <!-- comment -->\n"
         << "    <child1> text </child1>\n"
         << "    <child2> <child3 att=\"value\"/> </child2>\n"
         << "</tag>";
    file.close();

    XmlDocument doc(file.filename());
    CHECK(doc.file() == file.filename());

    XmlElement& root = doc.root();
    CHECK(root.tag() == "tag");
    CHECK(root.text().empty());
    CHECK(root.line() == 1);
    CHECK(root.document() == &doc);

    std::string att_str; int att_int;
    CHECK(root.hasAttribute("att"));
    CHECK(root.getAttribute("att", att_str) == "100");
    CHECK(root.getAttribute("att", att_int) == 100);
    CHECK_FALSE(root.hasAttribute("bla"));

    XmlElement::const_iterator it = root.begin();
    CHECK(it->tag() == "child1");
    CHECK(it->text() == " text ");
    CHECK(it->line() == 3);
    CHECK(it->begin() == it->end());

    it++;
    CHECK(it->tag() == "child2");
    CHECK(it->text().empty());
    CHECK(it->line() == 4);
    CHECK(it+1 == root.end());

    it = it->begin();
    CHECK(it->tag() == "child3");
    CHECK(it->text().empty());
    CHECK(it->line() == 4);
    CHECK(it->getAttribute("att", att_str) == "value");

    // Parsing errors
    CHECK_THROWS_AS(XmlElement("a <tag/>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag"), Error);
    CHECK_THROWS_AS(XmlElement("<tag>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag att=\" />"), Error);
    CHECK_THROWS_AS(XmlElement("<tag att=\"\" >"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> text <tag>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> text </ta>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> text </tag"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> <child> </tag>"), Error);

}

