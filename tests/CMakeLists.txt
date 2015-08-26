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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(set_state test_set_state.cpp)
target_link_libraries(set_state mutation++ ${Boost_LIBRARIES})

add_executable(stefan_maxwell test_stefan_maxwell.cpp)
target_link_libraries(stefan_maxwell mutation++ ${Boost_LIBRARIES})

add_executable(transfer_source test_transfer_source.cpp)
target_link_libraries(transfer_source mutation++ ${Boost_LIBRARIES})

add_executable(wdot test_wdot.cpp)
target_link_libraries(wdot mutation++ ${Boost_LIBRARIES})

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
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
            COMMAND ${test} ${mixture} )
    endforeach(mixture)
endfunction()

# Make sure mixtures load
test_all_mixtures(checkmix)

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

# Fortran examples
if (BUILD_FORTRAN_WRAPPER)
    test_example(wrapper_test fortran)
endif()

# Transport property tests
test_all_mixtures(stefan_maxwell)

# Kinetics tests
test_all_mixtures(wdot)

# StateModel tests
test_all_mixtures(set_state)

# TransferModel tests
test_all_mixtures(transfer_source)

#####################################
# Setup comparison tests
#####################################
add_executable(compare comparison.cpp)
target_link_libraries(compare mutation++)

file(GLOB files "comparisons/*.dat")
foreach(file ${files})
    get_filename_component(name ${file} NAME_WE)
    add_test(
        NAME compare_${name} 
        COMMAND compare ${file} 
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data )
endforeach()

