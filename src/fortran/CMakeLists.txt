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
cmake_policy(SET CMP0022 NEW)

include(GNUInstallDirs)

add_library(mutation++_fortran
    SHARED 
        cwrapper.cpp 
        cwrapper_interface.f90
)
add_library(mutation++::fortran_interface ALIAS mutation++_fortran)

set_target_properties(mutation++_fortran 
    PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(mutation++_fortran
    PRIVATE
        mutation++
)
install(TARGETS mutation++_fortran EXPORT mutation++Targets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/mutationpp.mod 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

