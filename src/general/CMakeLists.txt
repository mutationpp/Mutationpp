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

add_sources(mutation++
    Mixture.cpp
    MixtureOptions.cpp
)

add_headers(mutation++
    GlobalOptions.h 
    mutation++.h 
    Mixture.h 
    MixtureOptions.h 
    Constants.h 
    Errors.h
)

# Build and install executables
set(mutation++_exes checkmix mppequil bprime)
foreach (exe ${mutation++_exes})
    add_executable(${exe} ${exe}.cpp)
    target_link_libraries(${exe} 
        PRIVATE
            mutation++
    )
    install(
        TARGETS ${exe}
        EXPORT mutation++Targets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
)
endforeach(exe ${mutation++_exes})




