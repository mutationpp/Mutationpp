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

# Change this path to your Mutation++ install path
MUTATIONPP_INSTALL = ../../install

CXX_FLAGS = -I $(MUTATIONPP_INSTALL)/include \
			-I $(MUTATIONPP_INSTALL)/include/mutation++ \
 			-O3 -g

# use .dylib for Mac, .so for Unix
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CXX      = g++
	CXX_LIBS = $(MUTATIONPP_INSTALL)/lib/libmutation++.so
endif
ifeq ($(UNAME_S),Darwin)
	CXX      = c++
	CXX_LIBS = $(MUTATIONPP_INSTALL)/lib/libmutation++.dylib
endif

EXAMPLES  = \
	equilibrium_air \
	O2_dissociation \
	air_diffusion_comparison

all : $(EXAMPLES)

equilibrium_air : equilibrium_air.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

O2_dissociation : O2_dissociation.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

air_diffusion_comparison : air_diffusion_comparison.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS) $<

clean :
	$(RM) *.o $(EXAMPLES)

