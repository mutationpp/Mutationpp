#include "catch.hpp"
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
    SECTION(filename) {
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
            INFO("Line:" << i + 4)

            // Set the state of the mixture
            Eigen::Map<Eigen::VectorXd>(v1,n1) = data.row(i).segment( 0,n1);
            Eigen::Map<Eigen::VectorXd>(v2,n2) = data.row(i).segment(n1,n2);
            mix.setState(v1, v2, var_set);

            // Now compute the function and compare with expected result
            CHECK(function->compare(mix, data.row(i).tail(result_size)));
        }

        delete function;
    }
}

TEST_CASE("Comparison tests", "[comparisons]")
{
    comparisonTest("lambdae_air11_3rd_order.dat");
    comparisonTest("lambdae_air5.dat");
    comparisonTest("lambdah_air11_CG.dat");
    comparisonTest("lambdah_air11_LDLT.dat");
    comparisonTest("lambdah_air11_Wilke.dat");
    comparisonTest("lambdah_air5_CG.dat");
    comparisonTest("lambdah_air5_LDLT.dat");
    comparisonTest("lambdah_air5_Wilke.dat");
    comparisonTest("lambdaint_air11_NASA-9.dat");
    comparisonTest("lambdaint_air11_RRHO.dat");
    comparisonTest("lambdaint_air5_NASA-9.dat");
    comparisonTest("lambdaint_air5_RRHO.dat");
    comparisonTest("sigma_air11_2nd_order.dat");
    comparisonTest("sigma_air5.dat");
    comparisonTest("viscosity_air11_CG.dat");
    comparisonTest("viscosity_air11_Gupta-Yos.dat");
    comparisonTest("viscosity_air11_LDLT.dat");
    comparisonTest("viscosity_air11_Wilke.dat");
    comparisonTest("viscosity_air5_CG.dat");
    comparisonTest("viscosity_air5_Gupta-Yos.dat");
    comparisonTest("viscosity_air5_LDLT.dat");
    comparisonTest("viscosity_air5_Wilke.dat");
}
