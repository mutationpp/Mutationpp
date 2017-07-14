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

# Require Boost.Test package 
find_package( Boost COMPONENTS unit_test_framework filesystem system REQUIRED)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
configure_file(Configuration.h.in Configuration.h)

add_custom_target(
    generate_comparison_tests
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/generate_comparison_tests.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set(test_sources
    comparison_tests.cpp
    run_tests.cpp
    test_mixtures.cpp
)

add_executable(run_tests ${test_sources})
add_dependencies(run_tests generate_comparison_tests)
target_link_libraries(run_tests mutation++)

add_executable(diffusion_matrix test_diffusion_matrix.cpp)
target_link_libraries(diffusion_matrix mutation++)

add_executable(dXidT test_dXidT.cpp)
target_link_libraries(dXidT mutation++)

add_executable(energies test_energies.cpp)
target_link_libraries(energies mutation++)

add_executable(set_state test_set_state.cpp)
target_link_libraries(set_state mutation++)

add_executable(stefan_maxwell test_stefan_maxwell.cpp)
target_link_libraries(stefan_maxwell mutation++)

add_executable(transfer_source test_transfer_source.cpp)
target_link_libraries(transfer_source mutation++)

add_executable(wdot test_wdot.cpp)
target_link_libraries(wdot mutation++)

add_executable(thermal_diff_ratios test_thermal_diff_ratios.cpp)
target_link_libraries(thermal_diff_ratios mutation++)

#####################################
# TESTS
#####################################
set(mixtures_to_test
    air5_RRHO_ChemNonEq1T
    air5_RRHO_ChemNonEqTTv
    air5_NASA-7_ChemNonEq1T
    air5_NASA-9_ChemNonEq1T
    air11_RRHO_ChemNonEq1T
    air11_RRHO_ChemNonEqTTv
    air11_NASA-7_ChemNonEq1T
    air11_NASA-9_ChemNonEq1T
    argon_CR_ChemNonEq1T
    argon_CR_ChemNonEqTTv)
    
function(test_all_mixtures test)
    foreach(mixture ${mixtures_to_test})
        add_test(
            NAME ${test}_${mixture}
            COMMAND ${test} ${mixture}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data )
    endforeach(mixture)
endfunction()

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

# Equilibrium tests
test_all_mixtures(dXidT)

# Transport property tests
test_all_mixtures(diffusion_matrix)
test_all_mixtures(stefan_maxwell)
test_all_mixtures(thermal_diff_ratios)

# Kinetics tests
test_all_mixtures(wdot)

# StateModel tests
test_all_mixtures(energies)
test_all_mixtures(set_state)

# TransferModel tests
test_all_mixtures(transfer_source)

add_executable(update_comparison update_comparison.cpp)
target_link_libraries(update_comparison mutation++)


