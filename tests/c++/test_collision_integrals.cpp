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
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Catch;
using namespace Eigen;

/**
 * Helper function for loading the Q11, Q12, Q13, and Q22 collision integrals
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
 * Tests the constant collision integral type.
 */
TEST_CASE("Test constant collision integrals", "[transport]")
{
    // Generate a fake collision integral database
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"constant\" value=\"10\" />"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));

    // Make sure that the warning would be displayed
    SharedPtr<CollisionIntegral> Q11 = pair.get("Q11");
    REQUIRE(Q11->loaded());

    // Check that indeed the value is given
    CHECK(Q11->compute(1000.0) == 10.0);
    CHECK(Q11->compute(2000.0) == 10.0);
}

/**
 * Tests the exp-poly integral type.
 */
TEST_CASE("Test exp-poly collision integrals", "[transport]")
{
    // Generate a fake collision integral database (all the same)
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"exp-poly\" units=\"K,Å-Å\">"
         << "            -5.0e-3 1.0e-1 -9.8e-1 6.45e+0 </Q11>"
         << "        <Q12 type=\"exp-poly\" units=\"K,Å-Å\">"
         << "            -5.0e-3 1.0e-1 -9.8e-1 6.45e+0 </Q12>"
         << "        <Q13 type=\"exp-poly\" units=\"K,Å-Å\">"
         << "            -5.0e-3 1.0e-1 -9.8e-1 6.45e+0 </Q13>"
         << "        <Q22 type=\"exp-poly\" units=\"K,Å-Å\">"
         << "            -5.0e-3 1.0e-1 -9.8e-1 6.45e+0 </Q22>"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));

    // Temperatures used for reference values
    VectorXd T(5); T << 300.0, 500.0, 1000.0, 2000.0, 5000.0;

    // Reference data
    // Columns (Q11, Q12, Q13, Q22 A^2), Rows (300, 500, 1000, 2000, 5000 K)
    MatrixXd ref(5,4); ref <<
        24.1865, 24.1865, 24.1865, 24.1865,
        20.527,  20.527,  20.527,  20.527,
        16.511,  16.511,  16.511,  16.511,
        13.2345, 13.2345, 13.2345, 13.2345,
        9.6612,  9.6612,  9.6612,  9.6612;

    checkHeavyCollisionIntegrals(pair, T, ref);
}

/**
 * Tests the "from A*,B*,C*" integral types.
 */
TEST_CASE("Test from A*, B*, C* collision integrals", "[transport]")
{
    // Generate a fake collision integral database (all the same)
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"constant\" value=\"1e-20\" />"
         << "        <Q12 type=\"from C*\" />"
         << "        <Q13 type=\"from B*\" />"
         << "        <Q22 type=\"from A*\" />"
         << "        <Ast type=\"table\" units=\"m-m\" >"
         << "            500 5000, 1.0 20.0 </Ast>"
         << "        <Bst type=\"table\" units=\"m-m\" >"
         << "            500 5000, 5.0 15.0 </Bst>"
         << "        <Cst type=\"table\" units=\"m-m\" >"
         << "            500 5000, 3.0 22.0 </Cst>"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));
    SharedPtr<CollisionIntegral> Ast = pair.get("Ast");
    SharedPtr<CollisionIntegral> Bst = pair.get("Bst");
    SharedPtr<CollisionIntegral> Cst = pair.get("Cst");

    // Temperatures used for reference values
    VectorXd T(5); T << 500.0, 1000.0, 2000.0, 3000.0, 5000.0;

    // Reference data
    // Columns (Q11, Q12, Q13, Q22 A^2)
    MatrixXd ref(5,4);
    for (int i = 0; i < 5; ++i) {
        ref(i,0) = 1.0;
        ref(i,1) = Cst->compute(T(i));
        ref(i,2) = 1.25*Cst->compute(T(i)) - 0.25*Bst->compute(T(i));
        ref(i,3) = Ast->compute(T(i));
    }

    checkHeavyCollisionIntegrals(pair, T, ref);
}

/**
 * Tests the ratio integral type.
 */
TEST_CASE("Test ratio collision integrals", "[transport]")
{
    // Generate a fake collision integral database (all the same)
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"table\" units=\"K,Å-Å\">"
         << "            500 5000, 1.0 20.0 </Q11>"
         << "        <Q12 type=\"ratio\" ratio=\"1.0\" integral=\"Q11\" />"
         << "        <Q13 type=\"ratio\" ratio=\"2.0\" integral=\"Q11\" />"
         << "        <Q22 type=\"ratio\" ratio=\"3.0\" integral=\"Q11\" />"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));

    // Temperatures used for reference values
    VectorXd T(5); T << 500.0, 1000.0, 2000.0, 3000.0, 5000.0;

    // Reference data
    // Columns (Q11, Q12, Q13, Q22 A^2), Rows (500, 1000, 2000, 3000, 5000 K)
    MatrixXd ref(5,4); ref <<
        1.00000, 1.00000, 2.00000, 3.00000,
        3.11111, 3.11111, 6.22222, 9.33333,
        7.33333, 7.33333, 14.6667, 22.0000,
        11.5556, 11.5556, 23.1111, 34.6667,
        20.0000, 20.0000, 40.0000, 60.0000;

    checkHeavyCollisionIntegrals(pair, T, ref);
}

/**
 * Tests the Murphy integral type.
 */
TEST_CASE("Test Murphy collision integrals", "[transport]")
{
    // Generate a fake collision integral database (all the same)
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"table\" units=\"K,Å-Å\">"
         << "            500 5000, 1.0 20.0 </Q11>"
         << "        <Q12 type=\"table\" units=\"K,Å-Å\">"
         << "            500 5000, 2.0 15.0 </Q12>"
         << "        <Q13 type=\"Murphy\">"
         << "            <Q1 type=\"ratio\" ratio=\"1\" integral=\"Q11\" />"
         << "            <Q2 type=\"ratio\" ratio=\"1\" integral=\"Q11\" /> </Q13>"
         << "        <Q22 type=\"Murphy\">"
         << "            <Q1 type=\"ratio\" ratio=\"1\" integral=\"Q11\" />"
         << "            <Q2 type=\"ratio\" ratio=\"1\" integral=\"Q12\" /> </Q22>"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));
    SharedPtr<CollisionIntegral> Q11 = pair.get("Q11");
    SharedPtr<CollisionIntegral> Q12 = pair.get("Q12");

    // Temperatures used for reference values
    VectorXd T(5); T << 500.0, 1000.0, 2000.0, 3000.0, 5000.0;

    // Reference data
    // Columns (Q11, Q12, Q13, Q22 A^2), Rows (500, 1000, 2000, 3000, 5000 K)
    MatrixXd ref(5,4);
    for (int i = 0; i < 5; ++i) {
        ref(i,0) = Q11->compute(T(i));
        ref(i,1) = Q12->compute(T(i));
        ref(i,2) = std::sqrt(2.0) * Q11->compute(T(i)); // sqrt(Q11^2 + Q11^2)
        ref(i,3) = std::sqrt(ref(i,0)*ref(i,0) + ref(i,1)*ref(i,1));
    }
    ref *= 1e20;

    checkHeavyCollisionIntegrals(pair, T, ref);
}

/**
 * Tests the Pirani potential.
 */
TEST_CASE("Test Pirani potential collision integrals", "[transport]")
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

    // Create necessary species to define the pairs
    Species N2("N2"), O2("O2"), N2p("N2+");

    // Temperatures used for reference values
    VectorXd T(5); T << 300.0, 500.0, 1000.0, 2000.0, 5000.0;

    SECTION("Neutral-Neutral, explicit data") {
        // O2-O2 is defined with explicit beta, re, and phi0 values
        CollisionPair pair(O2, O2, &(doc.root()));

        // Reference data
        // Columns (Q11, Q12, Q13, Q22 A^2), Rows (300, 500, 1000, 2000, 5000 K)
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
        // Columns (Q11, Q12, Q13, Q22 A^2), Rows (300, 500, 1000, 2000, 5000 K)
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
        // Columns (Q11, Q12, Q13, Q22 A^2), Rows (300, 500, 1000, 2000, 5000 K)
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
        // Columns (Q11, Q12, Q13, Q22 A^2), Rows (300, 500, 1000, 2000, 5000 K)
        MatrixXd ref(5,4); ref <<
            97.82063,    78.77619,    66.65381,   102.23491,
            71.71063,    56.70362,    47.48123,    77.55206,
            46.45066,    36.37620,    30.96292,    52.13261,
            31.30984,    25.71801,    23.19851,    35.20541,
            22.16644,    20.17746,    19.16061,    24.56247;

        checkHeavyCollisionIntegrals(pair, T, ref);
    }
}

/**
 * Tests the table collision integral type.
 */
TEST_CASE("Test table collision integrals", "[transport]")
{
    // Generate a fake collision integral database which tests clipping and
    // interpolation
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"table\" units=\"K,Å-Å\">" // defaults to linear
         << "            1000  2000  4000  5000  6000,"
         << "            0.12  0.19  0.43  0.55  0.66 </Q11>"
         << "        <Q12 type=\"table\" units=\"K,Å-Å\" interpolator=\"MonotoneCubic\">"
         << "            1000  2000  4000  5000  6000,"
         << "            1.12  1.19  1.43  1.55  1.66 </Q12>"
         << "        <Q13 type=\"table\" units=\"K,Å-Å\" interpolator=\"Chebyshev\">"
         << "            1000  2000  4000  5000  6000,"
         << "            2.12  2.19  2.43  2.55  2.66 </Q13>"
         << "        <Q22 type=\"table\" units=\"K,Å-Å\" clip=\"no\">"
         << "            1000  2000  4000  5000  6000,"
         << "            3.12  3.19  3.43  3.55  3.66 </Q22>"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));

    // Temperatures used for reference values
    VectorXd T(4); T << 500.0, 1000.0, 3000.0, 5000.0;

    // Reference data
    // Columns (Q11, Q12, Q13, Q22 A^2), Rows (500, 1000, 3000, 5000 K)
    MatrixXd ref(4,4); ref <<
        0.12,  1.12,    2.12,    3.085,
        0.12,  1.12,    2.12,    3.12,
        0.31,  1.3015,  2.3053,  3.31,
        0.55,  1.55,    2.5485,  3.55;

    checkHeavyCollisionIntegrals(pair, T, ref);
}

/**
 * Tests the warning collision integral type.
 */
TEST_CASE("Test warning collision integrals", "[transport]")
{
    // Generate a fake collision integral database
    TemporaryFile file;
    file << "<collisions>"
         << "    <pair s1=\"N2\" s2=\"N2\">"
         << "        <Q11 type=\"warning\" value=\"10\" />"
         << "    </pair>"
         << "</collisions>";
    file.close();

    // Load the temporary collision database
    XmlDocument doc(file.filename());

    // Create necessary species to define the pairs
    Species N2("N2");

    // Create the collision pair
    CollisionPair pair(N2, N2, &(doc.root()));

    // Capture the standard output
    std::ostringstream oss;
    std::streambuf* cout_streambuf = std::cout.rdbuf();
    std::cout.rdbuf(oss.rdbuf());

    // Make sure that the warning would be displayed
    SharedPtr<CollisionIntegral> Q11 = pair.get("Q11");
    REQUIRE(Q11->loaded());

    // Restore standard output
    std::cout.rdbuf(cout_streambuf);

    // Check what was written
    CHECK(oss.str() ==
        "Warning: missing collision integral Q11_(N2,N2).  "
        "Using a constant value of 10.\n");

    // Check that indeed the value is given
    CHECK(Q11->compute(1000.0) == 10.0);
}
