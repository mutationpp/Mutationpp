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
    JacobianManager.cpp
    Kinetics.cpp
    RateLaws.cpp
    RateManager.cpp
    Reaction.cpp
    ReactionType.cpp
    StoichiometryManager.cpp
)

install(FILES JacobianManager.h DESTINATION include/mutation++)
install(FILES Kinetics.h DESTINATION include/mutation++)
install(FILES RateLaws.h DESTINATION include/mutation++)
install(FILES RateLawGroup.h DESTINATION include/mutation++)
install(FILES RateManager.h DESTINATION include/mutation++)
install(FILES Reaction.h DESTINATION include/mutation++)
install(FILES ReactionType.h DESTINATION include/mutation++)
install(FILES StoichiometryManager.h DESTINATION include/mutation++)
install(FILES ThirdBodyManager.h DESTINATION include/mutation++)

