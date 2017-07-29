#
# Copyright 2014-2015 von Karman Institute for Fluid Dynamics (VKI)
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

set(test_sources
    run_tests.cpp
    test_comparisons.cpp
    test_diffusion_matrix.cpp
    test_dXidT.cpp
    test_energies.cpp
    test_mixtures.cpp
    test_set_state.cpp
    test_species.cpp
    test_stefan_maxwell.cpp
    test_thermal_diff_ratios.cpp
    test_thermodb.cpp
    test_transfer_source.cpp
    test_wdot.cpp
)

add_executable(run_tests ${test_sources})
target_link_libraries(run_tests mutation++)

add_test(NAME Catch_tests COMMAND run_tests)

# Make sure examples compile and run
function(test_example name subdir)    
    add_test(
        NAME example_${name}_compile
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/${subdir}
        COMMAND make -B ${name} )
    add_test(
        NAME example_${name}_run
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/${subdir}
        COMMAND ${name} )
endfunction()

# C++ examples
test_example(equilibrium_air c++)
test_example(O2_dissociation c++)
test_example(air_diffusion_comparison c++)

# Fortran examples
if (BUILD_FORTRAN_WRAPPER)
    test_example(wrapper_test fortran)
endif()

add_executable(update_comparison update_comparison.cpp)
target_link_libraries(update_comparison mutation++)


