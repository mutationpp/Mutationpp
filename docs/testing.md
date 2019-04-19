# Mutation++ Testing Framework

**Contents** <br>
- [1. Overview](#1-overview-of-the-mutation-testing-framework)
  - [1.1 Catch](#11-catch)
  - [1.2 CTest](#12-ctest)
  - [1.3 Continuous integration with TravisCI](#13-continuous-integration-with-travis-ci)
- [2. Running the tests](#2-running-the-tests)
  - [2.1 Compile tests](#21-compile-tests)
  - [2.2 Using CTest](#22-using-ctest)
  - [2.3 Running specific tests](#23-running-specific-tests)
- [3. Writing new tests](#3-writing-new-tests)
  - [3.1 Comparison tests](#31-comparison-tests)
  - [3.2 Unit tests](#32-unit-tests)
  - [3.3 Testing example codes](#33-testing-example-codes)
- [4. Test coverage](#4-test-coverage)

# 1. Overview
[(TOC)](#table-of-contents)

Testing is an integral part of the software development process.  This is even more true for collaborative open source projects such as Mutation++ which uses a [continuous delivery system](https://en.wikipedia.org/wiki/Continuous_delivery), in which code is routinely updated, added, and deprecated.  This short guide will walk you through the testing framework used in Mutation++. Whether you are a developer or simply want to use the library as a [black box](https://en.wikipedia.org/wiki/Black_box), this guide will help you understand how tests are written and run.  For developers wishing to merge their code into main repository, you must follow these guidelines before any merge request is accepted.

The testing framework is composed of 3 main parts: the [Catch C++ unit testing library](https://github.com/philsquared/Catch), the CMake testing framework [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest), and finally continuous integration using [Travis CI](https://travis-ci.org/).  How each of these pieces fit into the Mutation++ testing framework is presented in the following subsections.

## 1.1 Catch
Catch is a single header C++ testing framework, designed to simplify writing C++ unit tests.  As a single header library, it provides an advantage over larger libraries in that it does not place any additional dependencies on Mutation++.  Catch represents the front line in the Mutation++ testing framework.  All unit and regression tests are implemented using the Catch macros.  Writing new tests using Catch is demonstrated in [section 3.2](#32-unit-tests).

## 1.2 CTest
CTest is the testing framework tool provided by the [CMake build environment](https://cmake.org).  CMake is an open-source cross-platform family of tools for building and deploying software and is used to generate the appropriate build files when [compiling and installing](installation-guide) Mutation++.  CTest provides an additional capability to run a range of tests, configured for the particular system that Mutation++ is installed on.  It is used to combine the unit and regression tests with other types of tests in Mutation++ into a single, coherent framework.

## 1.3 Continuous integration with Travis CI
Finally, Travis CI provides the ability to perform continuous integration (CI)..  This is the last piece of the Mutation++ testing puzzle.  Whenever a merge request is submitted, full test suite is run on a number of operating systems using multiple build configurations.  If the CI fails at any point during this process, a merge request cannot be accepted, thus preserving the stability of the master branch.  This is a crucial component of continuous delivery.  Of course, if something is not tested, then this process may still fail to catch bugs or other unintended consequences as they arise.  It is therefore up to developers and merge request reviewers to ensure that new code is properly tested before merging.

# 2. Running the tests 
[(TOC)](#table-of-contents)

As mentioned above, the testing framework is always run before accepting any changes into the library.  However, it is often useful to run the tests locally on your own machine.  Even if you are just a user, the best way to ensure that the library has been installed correctly is to run the tests after installation.  For developers, we recommend that you use a [test driven development](https://en.wikipedia.org/wiki/Test-driven_development) approach when writing new code, and thus it is necessary to run the tests several times as you implement your changes.

## 2.1 Compile tests
The first step in running the testing framework is to compile it.  By default, the testing framework is not installed, so this must be changed in the CMake configuration step by turning `ENABLE_TESTING` to `ON`.  Once this is done, recompile and install
```
make install
```
You should now see additional files being compiled related to the testing.  

## 2.2 Using CTest
Once compilation complete, run the tests using
```
ctest
```
in the `build` directory.  If all tests are successful, then the final output should appear as the following:

```
    Start 1: Catch_tests
1/7 Test #1: Catch_tests ................................   Passed   10.60 sec
    Start 2: example_equilibrium_air_compile
2/7 Test #2: example_equilibrium_air_compile ............   Passed   13.74 sec
    Start 3: example_equilibrium_air_run
3/7 Test #3: example_equilibrium_air_run ................   Passed    0.06 sec
    Start 4: example_O2_dissociation_compile
4/7 Test #4: example_O2_dissociation_compile ............   Passed    1.87 sec
    Start 5: example_O2_dissociation_run
5/7 Test #5: example_O2_dissociation_run ................   Passed    1.79 sec
    Start 6: example_air_diffusion_comparison_compile
6/7 Test #6: example_air_diffusion_comparison_compile ...   Passed    2.43 sec
    Start 7: example_air_diffusion_comparison_run
7/7 Test #7: example_air_diffusion_comparison_run .......   Passed    0.05 sec

100% tests passed, 0 tests failed out of 7

Total Test time (real) =  30.65 sec
```
Note that they test may take several minutes to run, depending on the build configuration.  Testing in `Debug` mode should be expected to be much slower than in `Release` mode.  The first test in the output above is actually all of the unit and regression tests written using Catch.  The remaining tests ensure that the example programs provided with Mutation++ compile and execute correctly.  This ensures that example programs are updated in the event of [API](https://en.wikipedia.org/wiki/Application_programming_interface) changes.  See the next section for running individual unit tests.  Finally, when a test fails, it is often useful to run ctest in verbose mode.  This is done using the `-V` option.  Use the `-h` option for a more detailed description of CTest.

## 2.3 Running specific tests
It is often useful to run specific tests only.  For example, if one test fails, you may only want to run that test again in verbose mode, in order to reduce all the output.  It is also useful during development to only test things you are changing instead of running all the tests over and over.  All the unit and regression tests are packaged in a single executable.  After building the tests, this executable can be run using
```
./tests/run_tests
```
from the `build` directory.  Note that this is all that is done when CTest runs the first test "Catch_tests" in the output above.  In the event that all tests are successful, the output will look something like
```
===============================================================================
All tests passed (36894 assertions in 13 test cases)
```
which shows that all tests passed, which comprised of 36894 individual assertions in 13 test cases.

As will be shown in the next section on writing new tests, Catch provides a very useful framework for naming individual test cases as well as tagging each test case with a set of key words which may be used to group selections of tests.  To see the list of test cases and tags, use the `-l` option with the `run_tests` command.  At the time of writing this guide, the output of this option was
```
All available test cases:
  Comparison tests
      [comparisons]
  DiffusionMatrix yields diffusion fluxes which sum to zero
      [DiffusionMatrix][transport]
  Equilibrium mole fractions derivatives sum to zero
      [equilibrium][thermodynamics]
  Sum of species energies equals mixture energies
      [thermodynamics]
  Mixture options
      [loading][mixtures]
  Loading test mixtures
      [loading][mixtures]
  setState() converts rho*Em to Tm and vice a versa
      [thermodynamics]
  stefanMaxwell() yields diffusion fluxes which sum to zero
      [transport]
  stefanMaxwell() yields zero net conduction current
      [transport]
  Thermal diffusion ratios sum to zero
      [transport]
  Energy transfer source terms are zero in equilibrium
      [transfer]
  Species production rates sum to zero
      [kinetics]
  Species production rates are zero in equilibrium
      [kinetics]
13 test cases
```
In order to run an individual test, simply run `run_tests` providing the test name as an option.
```
./tests/run_tests "Loading test mixtures"
===============================================================================
All tests passed (80 assertions in 1 test case)
```
Alternatively, tests associated with a given tag can be run together.
```
./tests/run_tests [transport]
===============================================================================
All tests passed (3400 assertions in 4 test cases)
```
To see the other options available, use the `-h` or `--help` options.

# 3. Writing new tests
[(TOC)](#table-of-contents)

This section explains how to add new tests to the testing framework.  Ideally, tests should be added whenever new features are being added. 

## 3.1 Comparison tests
Comparison tests are probably the easiest way to add tests to Mutation++.  Essentially, comparison tests provide a simple form of [regression testing](https://en.wikipedia.org/wiki/Regression_testing), tailored specifically to the set state / get property paradigm of Mutation++.  In short, a comparison test consists of the following pseudo-code:
```pseudo
Load a given mixture
for all comparisons:
    Set the state of the mixture using the given input
    Compute the comparison function result
    If abs(result - expected) > tolerance:
       The test fails
The test succeeds
```

In order to simplify tests of this kind, a small framework has been setup in Mutation++.
The following steps describe how to add a new comparison test to this framework.

### Add a new comparison file to `tests/comparisons`
Each file in `tests/comparisons` specifies a single comparison test.  For example, the first few lines of `tests/comparisons/lambdah_air5_LDLT.dat` are as follows:
```
air5_RRHO_ChemNonEq1T
1 5 1
heavy_thermal_conductivity 1.0e-10
  4.6114780801e-48  6.8557093657e-25  ...  5.0000000000e+02  2.9620534266e-02
  1.4508691285e-23  4.5794422313e-12  ...  1.0000000000e+03  4.8952825517e-02
  1.9701329711e-15  7.8153943326e-08  ...  1.5000000000e+03  6.3363445079e-02
  2.1825915686e-11  9.5655736490e-06  ...  2.0000000000e+03  7.7582314639e-02
```
The first 3 lines of the comparison test file represent the header, and tell the framework how to run the test.  Specifically, the first line indicates which mixture this test is to be run on (located in `tests/data/mixtures`).  The second line specifies three integer values:
  1. The variable set index used in the `setState()` function which tells the `Mixture` object what variables are being used to set the state.
  2. The size of the first vector passed to `setState()`.
  3. The size of the second vector passed to `setState()`.

The third line specifies the _comparison function_ and the tolerance for comparison.  Finally, the remaining lines in the comparison test file represent (state vector, expected result) pairs.  In the example above, each line indicates the 5 species densities, the temperature, and the expected heavy particle thermal conductivity at that state (for a total of 7 columns).  Thus, the loop in the pseudo-code above, represents a loop over each (state, result) line in the comparison test file.

### Add a new comparison function (if not already available)
The comparison function specified in the comparison test file described above, indicates which function is called to compare with the given data.  Comparison functions are added in the `tests/CompareFunc.h` file using the `ADD_FUNCTION` macro.  For example, the `heavy_thermal_conductivity` function, shown in the example, is defined in `tests/CompareFunc.h` as:
```c++
ADD_FUNCTION(heavy_thermal_conductivity, v(0) = mix.heavyThermalConductivity())
```
The macro takes two arguments as input.  The first argument specifies the name of the function.  This must correspond to the name used in the comparison test file.  The second argument specifies how to compute the result of the function.  Since comparison functions can compare whole vectors, the result is communicated as a vector.  In this case, only a single element is needed because the comparison function only compares a single value (the heavy particle thermal conductivity).  _The size of the result vector is determined automatically from the number of columns included in the comparison test file._  For example, a comparison function which compares the species thermal diffusion ratios could be defined as
```c++
ADD_FUNCTION(thermal_diffusion_ratios, mix.thermalDiffusionRatios(v.data()))
```

### Add the comparison test to `tests/test_comparisons.cpp`
Lastly, the framework must be told to run the new comparison test.  This is done by adding a new line to `tests/test_comparisons.cpp` which specifies the appropriate comparison file.  For example,
```c++
comparisonTest("lambdah_air5_LDLT.dat");
```

And that's it! After recompiling, you should be able to run the newly added comparison test.

## 3.2 Unit tests
- Writing unit tests with Catch
- Useful macros

```c++
#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include "catch.hpp"
#include <Eigen/Dense>

using namespace Mutation;
using namespace Catch;
using namespace Eigen;

TEST_CASE
(
    "Equilibrium mole fractions derivatives sum to zero",
    "[equilibrium][thermodynamics]"
)
{
    const double tol = std::numeric_limits<double>::epsilon();

    MIXTURE_LOOP
    (
        VectorXd dx(mix.nSpecies());

        EQUILIBRATE_LOOP
        (
            mix.dXidT(dx.data());
            INFO("dx = " << dx);
            REQUIRE(dx.sum() == Approx(0.0).margin(tol));
        )
    )
}
```

## 3.3 Testing example codes
If an example is included in the `examples` directory, then at a bare minimum, it should be able to compile and run without producing an error.  This is done by adding a new test to `tests/CMakeLists.txt` using the `test_example` function.  The `test_example` function takes the name of the example code (the file name without the extension) and the sub-directory of `examples` in which the test is placed.  For example,

```cmake
test_example(equilibrium_air c++)
```

will make sure that `examples/c++/equilibrium_air.cpp` compiles and runs successfully.  Each call to `test_example` in `tests/CMakeLists.txt` adds two tests to the CTest test list named `example_{name}_compile` and `example_{name}_run` where `{name}` corresponds to the name of the example.

# 4. Test Coverage
It is useful to see how much of the library the testing actually covers in terms of lines of code and functions.  For this purpose, you can build a special build type which generates a make rule which generates a code coverage analysis.

To try the code coverage, use the following command in the `build` directory.
```
cmake -DCMAKE_BUILD_TYPE=Coverage ..
make -jN test_coverage
```
where `N` is the number of processors to use for parallel compilation.  Only the `Coverage` build type will provide the `test_coverage` target.  Once built, you can view the results by opening `coverage/index.html` from the build directory in your favourite web browser.  
