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
#include <catch.hpp>

using namespace Mutation;
using namespace Mutation::Kinetics;
using namespace Mutation::Utilities::IO;
using namespace Catch;

/**
 * @brief Checks that the Arrhenius rate law works as expected.
 */
void checkArrheniusRate(
    double A, double n, double theta, int order, bool use_Ea = false)
{
    std::stringstream xml;
    xml <<  "<arrhenius A=\"" << A << "\" n = \"" << n << "\" ";
    if (use_Ea)
        xml << "Ea=\"" << theta * RU;
    else
        xml << "T=\"" << theta;
    xml << "\" />";

    Arrhenius rate(XmlElement(xml.str()), 1);

    CHECK(rate.A() == Approx(A));
    CHECK(rate.n() == Approx(n));
    CHECK(rate.T() == Approx(theta));

    for (int i = 0; i < 10; ++i) {
        double T = i * 1000.0 + 1000.0;
        double k = A * std::pow(T, n) * std::exp(-theta / T);

        CHECK( rate.getLnRate(std::log(T), 1.0/T) == Approx(std::log(k)) );

        CHECK(
            rate.derivative(k, std::log(T), 1.0/T) ==
            Approx(k * (n + theta / T) / T)
        );
    }
}

/**
 * Tests if the Arrhenius rate law works correctly.
 */
TEST_CASE("Arrhenius rate laws work correctly", "[kinetics]")
{
    // Check that the rate is calculated correctly
    for (int order = 1; order < 4; ++order) {
        checkArrheniusRate(5.0e15,  0.0, order, 75500.0);
        checkArrheniusRate(5.0e15,  0.0, order, 75500.0, true);
        checkArrheniusRate(5.0e15,  1.0, order, 75500.0);
        checkArrheniusRate(5.0e15,  1.0, order, 75500.0, true);
        checkArrheniusRate(5.0e15, -1.0, order, 75500.0);
        checkArrheniusRate(5.0e15, -1.0, order, 75500.0, true);
    }

    // Check that loading fails properly
    CHECK_THROWS( Arrhenius(XmlElement("<arrhenius />"), 1) );
    CHECK_THROWS( Arrhenius(XmlElement("<arrhenius A=\"1\"/>"), 1) );
    CHECK_THROWS( Arrhenius(XmlElement("<arrhenius A=\"1\" n=\"1\"/>"), 1) );
    CHECK_THROWS( Arrhenius(XmlElement("<arrhenius T=\"1\" n=\"1\"/>"), 1) );
    CHECK_THROWS( Arrhenius(XmlElement("<arrhenius Ea=\"1\" n=\"1\"/>"), 1) );
    CHECK_NOTHROW( Arrhenius(XmlElement("<arrhenius A=\"1\" T=\"1\"/>"), 1) );
    CHECK_NOTHROW( Arrhenius(XmlElement("<arrhenius A=\"1\" Ea=\"1\"/>"), 1) );
}
