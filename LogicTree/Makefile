# Change this path to your Mutation++ install path
MUTATIONPP_INSTALL = /Users/jbscoggi/src/mpp_branches/aurelie/install

CXX       = g++
CXX_FLAGS = -I$(MUTATIONPP_INSTALL)/include/mutation++ -O3 -g
# use .dylib for Mac, .so for Unix
CXX_LIBS  = $(MUTATIONPP_INSTALL)/lib/libmutation++.dylib

TREE  = \
	tree \
        ReactionType
	
all : ReactionType

ReactionType: ReactionType.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS) $<

clean :
	$(RM) *.o $(TREE)

