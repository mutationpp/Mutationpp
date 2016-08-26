/**
 * @file comparison.cpp
 *
 * @brief Generic tool for comparing results with a file.
 */

/*
 * Copyright 2014-2015 von Karman Institute for Fluid Dynamics (VKI)
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

#include "CompareFunc.h"
#include "MatrixHelper.h"
using namespace Mutation;

#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

/**
 * Driver program for comparing a function to the values stored in the given
 * file.
 */
int main(int argc, char* argv[])
{
    // Open the file
    ifstream file(argv[1]);

    // Load the mixture
    string str;
    getline(file, str);
    Mixture mix(str);

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
    bool matches = true;

    for (int i = 0; i < data.rows(); ++i) {
        // Set the state of the mixture
        Map<VectorXd>(v1,n1) = data.row(i).segment( 0,n1);
        Map<VectorXd>(v2,n2) = data.row(i).segment(n1,n2);
        mix.setState(v1, v2, var_set);

        // Now compute the function and compare with expected result
        matches &= function->compare(mix, data.row(i).tail(result_size));
    }

    // Return 0 on success
    return ((int) !matches);
}
