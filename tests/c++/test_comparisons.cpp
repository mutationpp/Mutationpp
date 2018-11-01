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

#include <catch.hpp>
#include "mutation++.h"
#include "Configuration.h"
#include "CompareFunc.h"
#include "MatrixHelper.h"

#include <string>
#include <Eigen/Dense>

using namespace Catch;
using namespace Mutation;


void comparisonTest(const std::string& filename)
{
    // Get the folder that is configured at compile time by CMake
    std::string comparisons_folder = COMPARISONS_FOLDER;

    // Open the file
    std::ifstream file((comparisons_folder + filename).c_str());

    // If we don't find the input file, then the test should fail
    REQUIRE(file.is_open());

    // Load the mixture
    GlobalOptions::workingDirectory(TEST_DATA_FOLDER);
    std::string str;
    getline(file, str);
    Mixture mix(str);

    // Load the state variable information
    getline(file, str);
    int var_set, n1, n2;
    std::stringstream s1(str);  s1 >> var_set >> n1 >> n2;

    // Load the function for comparison and tolerance
    getline(file, str);
    std::string fun; double tol;
    std::stringstream s2(str); s2 >> fun >> tol;

    CompareFunc* function =
        Utilities::Config::Factory<CompareFunc>::create(fun, tol);

    // Load the state variables and data to compare to
    Eigen::MatrixXd data = readMatrix(file);
    const int result_size = data.cols()-n1-n2;

    // Close the file
    file.close();

    // Loop over the states
    double v1 [n1];
    double v2 [n2];

    for (int i = 0; i < data.rows(); ++i) {
        INFO("Line:" << i + 4);

        // Set the state of the mixture
        Eigen::Map<Eigen::VectorXd>(v1,n1) = data.row(i).segment( 0,n1);
        Eigen::Map<Eigen::VectorXd>(v2,n2) = data.row(i).segment(n1,n2);
        mix.setState(v1, v2, var_set);

        // Now compute the function and compare with expected result
        try {
            bool passed =
                function->compare(mix, data.row(i).tail(result_size));
            CHECK(passed);
        } catch (Error& e) {
            INFO(e.what());
            CHECK(false);
        }
    }

    delete function;
}

TEST_CASE("Comparison tests", "[comparisons]")
{
    SECTION("lambdae_air11_3rd_order.dat") {
        comparisonTest("lambdae_air11_3rd_order.dat"); }
    SECTION("lambdae_air5.dat") {
        comparisonTest("lambdae_air5.dat"); }
    SECTION("lambdah_air11_CG.dat") {
        comparisonTest("lambdah_air11_CG.dat"); }
    SECTION("lambdah_air11_LDLT.dat") {
        comparisonTest("lambdah_air11_LDLT.dat"); }
    SECTION("lambdah_air11_Wilke.dat") {
        comparisonTest("lambdah_air11_Wilke.dat"); }
    SECTION("lambdah_air5_CG.dat") {
        comparisonTest("lambdah_air5_CG.dat"); }
    SECTION("lambdah_air5_LDLT.dat") {
        comparisonTest("lambdah_air5_LDLT.dat"); }
    SECTION("lambdah_air5_Wilke.dat") {
        comparisonTest("lambdah_air5_Wilke.dat"); }
    SECTION("lambdaint_air11_NASA.dat") {
        comparisonTest("lambdaint_air11_NASA-9.dat"); }
    SECTION("lambdaint_air11_RRHO.dat") {
        comparisonTest("lambdaint_air11_RRHO.dat"); }
    SECTION("lambdaint_air5_NASA.dat") {
        comparisonTest("lambdaint_air5_NASA-9.dat"); }
    SECTION("lambdaint_air5_RRHO.dat") {
        comparisonTest("lambdaint_air5_RRHO.dat"); }
    SECTION("sigma_air11_2nd_order.dat") {
        comparisonTest("sigma_air11_2nd_order.dat"); }
    SECTION("sigma_air5.dat") {
        comparisonTest("sigma_air5.dat"); }
    SECTION("viscosity_air11_CG.dat") {
        comparisonTest("viscosity_air11_CG.dat"); }
    SECTION("viscosity_air11_Gupta.dat") {
        comparisonTest("viscosity_air11_Gupta-Yos.dat"); }
    SECTION("viscosity_air11_LDLT.dat") {
        comparisonTest("viscosity_air11_LDLT.dat"); }
    SECTION("viscosity_air11_Wilke.dat") {
        comparisonTest("viscosity_air11_Wilke.dat"); }
    SECTION("viscosity_air5_CG.dat") {
        comparisonTest("viscosity_air5_CG.dat"); }
    SECTION("viscosity_air5_Gupta.dat") {
        comparisonTest("viscosity_air5_Gupta-Yos.dat"); }
    SECTION("viscosity_air5_LDLT.dat") {
        comparisonTest("viscosity_air5_LDLT.dat"); }
    SECTION("viscosity_air5_Wilke.dat") {
        comparisonTest("viscosity_air5_Wilke.dat"); }
}
