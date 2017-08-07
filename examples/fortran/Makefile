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

FC       = gfortran
FC_FLAGS = -J$(MUTATIONPP_INSTALL)/include/mutation++ -O3 -g -fdefault-real-8

# use .dylib for Mac, .so for Unix
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	FC_LIBS = $(MUTATIONPP_INSTALL)/lib/libmutation++_fortran.so
endif
ifeq ($(UNAME_S),Darwin)
	FC_LIBS = $(MUTATIONPP_INSTALL)/lib/libmutation++_fortran.dylib
endif

EXAMPLES  = \
	wrapper_test

all : $(EXAMPLES)

wrapper_test : wrapper_test.o
	$(FC) -o $@ $(FC_FLAGS) $^ $(FC_LIBS)

%.o : %.f90
	$(FC) -c $(FC_FLAGS) $<

clean :
	$(RM) *.o $(EXAMPLES)

