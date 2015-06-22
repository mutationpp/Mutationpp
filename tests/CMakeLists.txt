#
# Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

add_executable(gsi_catalysis test_gas_surface_interaction.cpp)
target_link_libraries(gsi_catalysis mutation++ ${Boost_LIBRARIES})

add_executable(catalysis_wall_solver test_catalysis_wall_solver.cpp)
target_link_libraries(catalysis_wall_solver mutation++ ${Boost_LIBRARIES})

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
    
set(wall_solver_mixture_to_test
    O2_gamma_001)
#    O2_gamma_005
#    O2_gamma_01
#    O2_gamma_05
#    O2_gamma_1)

function(test_all_mixtures test)
    foreach(mixture ${mixtures_to_test})
        add_test(
            NAME ${test}_${mixture}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
            COMMAND ${test} ${mixture} )
    endforeach(mixture)
endfunction()

function(test_gsi gsi_type)
       add_test(
           NAME gsi_${gsi_type}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/gsi
           COMMAND gsi_${gsi_type} )
   foreach(wall_solver_mixture ${wall_solver_mixture_to_test})
       add_test(
           NAME ${gsi_type}_wall_solver_${wall_solver_mixture_to_test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/gsi/wall_solver
           COMMAND ${gsi_type}_wall_solver ${wall_solver_mixture_to_test} )
   endforeach(wall_solver_mixture)
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

test_example(equilibrium_air g++)
test_example(O2_dissociation g++)
test_example(wrapper_test fortran)

# Transport property tests
test_all_mixtures(stefan_maxwell)

# Kinetics tests
test_all_mixtures(wdot)

# StateModel tests
test_all_mixtures(set_state)

# TransferModel tests
test_all_mixtures(transfer_source)

# Catalysis tests
test_gsi(catalysis)
