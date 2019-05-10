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

cmake_minimum_required(VERSION 2.8.2)

# Append the PROJECT_BINARY_DIR to CMAKE_PREFIX_PATH in order to be able to 
# find the configuration files at build time
list(APPEND CMAKE_PREFIX_PATH ${PROJECT_BINARY_DIR}/src)

# c++ examples
add_subdirectory(c++)

# Fortran examples
if (BUILD_FORTRAN_WRAPPER)
    add_subdirectory(fortran)
endif(BUILD_FORTRAN_WRAPPER)
