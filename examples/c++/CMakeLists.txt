# Macro to list all the subdirs of a dir
# Copyright 2014-2017 von Karman Institute for Fluid Dynamics (VKI)
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

cmake_minimum_required(VERSION 2.6)

include(ListSubDirs)
ListSubDirs(example_subdirs ${CMAKE_CURRENT_LIST_DIR})

# Add all examples targets and test them
foreach(subdir ${example_subdirs})
    add_subdirectory(${subdir})
    add_test(
        NAME example_${subdir}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${subdir}
        COMMAND ${subdir} )
endforeach()
