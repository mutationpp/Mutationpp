#
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
#

cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    ChemNonEqStateModel.cpp
    ChemNonEqTTvStateModel.cpp
    Composition.cpp
    EquilStateModel.cpp
    MultiPhaseEquilSolver.cpp
    Nasa7Polynomial.cpp
    Nasa9Polynomial.cpp
    Nasa7DB.cpp
    Nasa9DB.cpp
    ParticleRRHO.cpp
    RrhoDB.cpp
    Species.cpp
    SpeciesListDescriptor.cpp
    SpeciesNameFSM.cpp
    Thermodynamics.cpp
    ThermoDB.cpp
)

install(FILES Composition.h DESTINATION include/mutation++)
install(FILES MultiPhaseEquilSolver.h DESTINATION include/mutation++)
install(FILES ParticleRRHO.h DESTINATION include/mutation++)
install(FILES Species.h DESTINATION include/mutation++)
install(FILES SpeciesListDescriptor.h DESTINATION include/mutation++)
install(FILES SpeciesNameFSM.h DESTINATION include/mutation++)
install(FILES StateModel.h DESTINATION include/mutation++)
install(FILES ThermoDB.h DESTINATION include/mutation++)
install(FILES Thermodynamics.h DESTINATION include/mutation++)
