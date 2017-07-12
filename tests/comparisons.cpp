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

#include "catch.hpp"
#include <iostream>
#include <string>
#include "mutation++.h"
#include <Eigen/Dense>
#include "MatrixHelper.h"
#include "CompareFunc.h"
#include "Configuration.h"

using namespace std;
using namespace Catch;
using namespace Mutation;
using namespace Eigen;

/// Driver for the comparison tests.
void run_comparison(const char* filename)
{
    // Get the folder that is configured at compile time by CMake
    std::string comparisons_folder = COMPARISONS_FOLDER;

    // Open the file
    std::ifstream file((comparisons_folder + filename).c_str());

    // If we don't find the input file, then the test should fail
    REQUIRE(file.is_open());

    // Load the mixture
    std::string data_folder = TEST_DATA_FOLDER;
    std::string str;
    getline(file, str);
    Mixture mix(data_folder + str);

    // Load the state variable information
    getline(file, str);
    int var_set, n1, n2;
    stringstream s1(str);  s1 >> var_set >> n1 >> n2;

    // Load the function for comparison and tolerance
    getline(file, str);
    string fun; double tol;
    stringstream s2(str); s2 >> fun >> tol;

    CompareFunc* function =
        Utilities::Config::Factory<CompareFunc>::create(fun, tol);

    // Load the state variables and data to compare to
    MatrixXd data = readMatrix(file);
    const int result_size = data.cols()-n1-n2;

    // Close the file
    file.close();

    // Loop over the states
    double v1 [n1];
    double v2 [n2];

    for (int i = 0; i < data.rows(); ++i) {
        // Set the state of the mixture
        Map<VectorXd>(v1,n1) = data.row(i).segment( 0,n1);
        Map<VectorXd>(v2,n2) = data.row(i).segment(n1,n2);
        mix.setState(v1, v2, var_set);

        // Now compute the function and compare with expected result
        CHECK(function->compare(mix, data.row(i).tail(result_size)));
    }
}

//==============================================================================

TEST_CASE("Comparison: air11 3rd order electron thermal conductivity", "[transport][comparison]") {
    run_comparison("lambdae_air11_3rd_order.dat");
}
