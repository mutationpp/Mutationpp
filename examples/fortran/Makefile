# Change this path to your Mutation++ install path
MUTATIONPP_INSTALL = ../../install

FC       = gfortran
FC_FLAGS = -J$(MUTATIONPP_INSTALL)/include/mutation++ -O3 -g -fdefault-real-8
# use .dylib for Mac, .so for Unix
FC_LIBS  = $(MUTATIONPP_INSTALL)/lib/libmutation++_fortran.dylib

EXAMPLES  = \
	wrapper_test

all : $(EXAMPLES)

wrapper_test : wrapper_test.o
	$(FC) -o $@ $(FC_FLAGS) $^ $(FC_LIBS)

%.o : %.f90
	$(FC) -c $(FC_FLAGS) $<

clean :
	$(RM) *.o $(EXAMPLES)

