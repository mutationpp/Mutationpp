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
#include <catch.hpp>
#include <Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Catch;
using namespace Eigen;

/**
 * Helper function for loading the Q11, Q12, Q13, and Q14 collision integrals
 * for a given pair and comparing the computed integral values agains a given
 * set of reference data.
 */
void checkHeavyCollisionIntegrals(
    CollisionPair& pair, const VectorXd& T, const MatrixXd& ref)
{
    REQUIRE(T.size() == ref.rows());

    SharedPtr<CollisionIntegral> Q11 = pair.get("Q11"); REQUIRE(Q11->loaded());
    SharedPtr<CollisionIntegral> Q12 = pair.get("Q12"); REQUIRE(Q12->loaded());
    SharedPtr<CollisionIntegral> Q13 = pair.get("Q13"); REQUIRE(Q13->loaded());
    SharedPtr<CollisionIntegral> Q22 = pair.get("Q22"); REQUIRE(Q22->loaded());

    for (int i = 0; i < T.size(); ++i) {
        CHECK(Q11->compute(T(i))*1.0e20 == Approx(ref(i,0)));
        CHECK(Q12->compute(T(i))*1.0e20 == Approx(ref(i,1)));
        CHECK(Q13->compute(T(i))*1.0e20 == Approx(ref(i,2)));
        CHECK(Q22->compute(T(i))*1.0e20 == Approx(ref(i,3)));
    }
}

/**
 * Tests the Pirani potential.
 */
TEST_CASE
(
    "Pirani potential collision integrals",
    "[transport]"
)
{
    // Generate a fake collision integral database
    // From Laricchiuta2007
    // O2-N2:  beta = 8.1055, eps0 = 11.66, re = 3.805
    // O2-O2:  beta = 8.1375, eps0 = 11.97, re = 3.78
    // N2-N2+: beta = 8.0746, eps0 = 96.57, re = 3.454
    // O2-N2+: beta = 8.1055, eps0 = 90.61, re = 3.434
    TemporaryFile file;
    file << "<collisions>"
         << "    <dipole-polarizabilities units=\"Å-Å-Å\">"
         << "        <species name=\"N2\" value=\"1.74\" />"
         << "        <species name=\"N2+\" value=\"1.75\" />"
         << "        <species name=\"O2\" value=\"1.60\" />"
         << "    </dipole-polarizabilities>"
         << "    <effective-electrons>"
         << "        <species name=\"N2\" value=\"7.6\" />"
         << "        <species name=\"O2\" value=\"9.3333333333\" />"
         << "    </effective-electrons>"
         << "    <pair s1=\"O2\" s2=\"O2\">"
         << "        <Pirani-Data beta=\"8.1375\" re=\"3.78\" phi0=\"11.97\"/>"
         << "        <Q11 type=\"Pirani\"/>"
         << "        <Q12 type=\"Pirani\"/>"
         << "        <Q13 type=\"Pirani\"/>"
         << "        <Q22 type=\"Pirani\"/>"
         << "    </pair>"
         << "    <pair s1=\"N2\" s2=\"O2\">"
         << "        <Q11 type=\"Pirani\"/>"
         << "        <Q12 type=\"Pirani\"/>"
         << "        <Q13 type=\"Pirani\"/>"
         << "        <Q22 type=\"Pirani\"/>"
         << "    </pair>"
         << "    <pair s1=\"O2\" s2=\"N2+\">"
         << "        <Pirani-Data beta=\"8.1055\" re=\"3.434\" phi0=\"90.61\"/>"
         << "        <Q11 type=\"Pirani\"/>"
         << "        <Q12 type=\"Pirani\"/>"
         << "        <Q13 type=\"Pirani\"/>"
         << "        <Q22 type=\"Pirani\"/>"
         << "    </pair>"
         << "    <pair s1=\"N2+\" s2=\"N2\">"
         << "        <Q11 type=\"Pirani\"/>"
         << "        <Q12 type=\"Pirani\"/>"
         << "        <Q13 type=\"Pirani\"/>"
         << "        <Q22 type=\"Pirani\"/>"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Temperatures used for reference values
    VectorXd T(5); T << 300.0, 500.0, 1000.0, 2000.0, 5000.0;

    // Create necessary species to define the pairs
    Species N2("N2"), O2("O2"), N2p("N2+");

    SECTION("Neutral-Neutral, explicit data") {
        // O2-O2 is defined with explicit beta, re, and phi0 values
        CollisionPair pair(O2, O2, &(doc.root()));

        // Reference data
        // Columns (Q11, Q12, Q13, Q14 A^2), Rows (300, 500, 1000, 2000, 5000 K)
        MatrixXd ref(5,4); ref <<
            36.33002,    31.38263,    28.98695,    40.62648,
            30.66071,    27.80487,    26.38351,    33.99039,
            26.43153,    24.85450,    23.79778,    29.47545,
            23.60094,    22.24057,    21.22414,    26.79771,
            19.93869,    18.62013,    17.68011,    23.08201;

        checkHeavyCollisionIntegrals(pair, T, ref);
    }

    SECTION("Neutral-Neutral, implicit data") {
        // N2-O2 is defined with implicit data using the species polarizability
        // and effective electrons
        CollisionPair pair(N2, O2, &(doc.root()));

        // Reference data
        // Columns (Q11, Q12, Q13, Q14 A^2), Rows (300, 500, 1000, 2000, 5000 K)
        MatrixXd ref(5,4); ref <<
            36.37455,    31.49629,    29.14517,    40.66233,
            30.78823,    27.98169,    26.57213,    34.13051,
            26.60832,    25.03100,    23.96368,    29.69713,
            23.75599,    22.37806,    21.34904,    26.99934,
            20.03925,    18.70630,    17.75762,    23.21459;

        checkHeavyCollisionIntegrals(pair, T, ref);
    }

    SECTION("Ion-Neutral, explicit data") {
        // N2+-O2 is defined with explicit beta, re, and phi0 values
        CollisionPair pair(N2p, O2, &(doc.root()));

        // Reference data
        // Columns (Q11, Q12, Q13, Q14 A^2), Rows (300, 500, 1000, 2000, 5000 K)
        MatrixXd ref(5,4); ref <<
            93.35611,    75.03173,    63.39025,    97.95099,
            68.32777,    53.92690,    45.14383,    74.18166,
            44.29329,    34.74493,    29.69348,    49.79153,
            30.08704,    24.88870,    22.57074,    33.77278,
            21.61377,    19.76986,    18.79716,    23.96274;

        checkHeavyCollisionIntegrals(pair, T, ref);
    }

    SECTION("Ion-Neutral, implicit data") {
        // N2+-N2 is defined with implicit data using the species polarizability
        // and effective electrons
        CollisionPair pair(N2p, N2, &(doc.root()));

        // Reference data
        // Columns (Q11, Q12, Q13, Q14 A^2), Rows (300, 500, 1000, 2000, 5000 K)
        MatrixXd ref(5,4); ref <<
            97.82063,    78.77619,    66.65381,   102.23491,
            71.71063,    56.70362,    47.48123,    77.55206,
            46.45066,    36.37620,    30.96292,    52.13261,
            31.30984,    25.71801,    23.19851,    35.20541,
            22.16644,    20.17746,    19.16061,    24.56247;

        checkHeavyCollisionIntegrals(pair, T, ref);
    }
}

