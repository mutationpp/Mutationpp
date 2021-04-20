#
# Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
#
# This file is part of MUlticomponent Thermodynamic And Transport
# properties for IONized gases in C++ (Mutation++) software package.
#
# Mutation++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Mutation++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Mutation++.  If not, see
# <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
configure_file(Configuration.h.in Configuration.h)

file(GLOB test_sources "test_*.cpp")
add_executable(run_tests run_tests.cpp ${test_sources})
target_link_libraries(run_tests 
    PRIVATE
        mutation++
    PUBLIC
        Catch2::Catch2
)

# Add all Catch TEST_CASEs as individual CTest tests
include(CTest)
include(ParseAndAddCatchTests)
set(PARSE_CATCH_TESTS_VERBOSE ON)
ParseAndAddCatchTests(run_tests)

add_executable(update_comparison update_comparison.cpp)
target_link_libraries(update_comparison mutation++)

add_test(
    NAME bugfix_120
    COMMAND mppequil -T 300 -P 101325 -s 15,16 --species-list "O O2"
)
