# Change this path to your Mutation++ install path
MUTATIONPP_INSTALL = ../../install

CXX       = g++
CXX_FLAGS = -I$(MUTATIONPP_INSTALL)/include/mutation++ -O3 -g 
# use .dylib for Mac, .so for Unix
CXX_LIBS  = $(MUTATIONPP_INSTALL)/lib/libmutation++.so

EXAMPLES  = \
	equilibrium_air \
	O2_dissociation \
	Reaction_Rate_T_Tv_Te\
	rhoie_test

OUTPUT_FILES = \
	Reaction_Rates.dat

all : $(EXAMPLES)

equilibrium_air : equilibrium_air.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

O2_dissociation : O2_dissociation.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

Reaction_Rate_T_Tv_Te : Reaction_Rate_T_Tv_Te.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

rhoie_test : rhoie_test.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS) $<

clean :
	$(RM) *.o $(EXAMPLES) $(OUTPUT_FILES)

