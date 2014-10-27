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

add_executable(stefan_maxwell test_stefan_maxwell.cpp)
target_link_libraries(stefan_maxwell mutation++ ${Boost_LIBRARIES})
add_executable(wdot test_wdot.cpp)
target_link_libraries(wdot mutation++ ${Boost_LIBRARIES})

# Define the tests to be run
add_test(checkmix_air5 checkmix air5)
add_test(checkmix_air11 checkmix air11)
add_test(stefan-maxwell_air5 stefan_maxwell air5)
add_test(stefan-maxwell_air11 stefan_maxwell air11)
add_test(wdot_air5  wdot air5)
add_test(wdot_air11 wdot air11)

