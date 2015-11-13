#
# Copyright 2015 von Karman Institute for Fluid Dynamics (VKI)
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
 
# Require Boost.Python package 
find_package(Boost COMPONENTS python REQUIRED)
message(STATUS "Boost Python library is: ${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Require Python
find_package(PythonLibs 2.7 REQUIRED)
# Hack to force boost to use python2.7, needs to be fixed for production release
string (REGEX REPLACE "[0-9].[0-9]m" "2.7" PYTHON_LIBRARIES ${PYTHON_LIBRARIES}) 
string (REGEX REPLACE "[0-9].[0-9]m" "2.7" PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
message (STATUS "Python library is: ${PYTHON_LIBRARIES}")
message (STATUS "Python include dir is: ${PYTHON_INCLUDE_DIRS}")
link_directories(${PYTHON_LIBRARIES})
include_directories(${PYTHON_INCLUDE_DIRS})

# Build pyMutation
add_library(pyMutation SHARED
	numpy/dtype.cpp
	numpy/scalars.cpp
	numpy/ndarray.cpp
	numpy/matrix.cpp
	numpy/ufunc.cpp
	numpy/numpy.cpp
	pyMutation.cpp 
	export_MixtureOptions.cpp
	export_Mixture.cpp)

# Link mutation++, boost.python and python libraries
target_link_libraries(pyMutation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} mutation++)
