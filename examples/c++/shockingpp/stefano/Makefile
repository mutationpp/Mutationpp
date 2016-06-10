# CC = g++ 
# No optimization (to help Valgrind)
CC = g++ -g -O0
CXXFLAGS = 
LDFLAGS =


MUTATIONPP_INSTALL = /home/starlight/mutation++/mutation++/install
BOOSTINSTALL = /home/starlight/Downloads/boost_1_54_0

INCLUDE = -I$(BOOSTINSTALL)
INCLUDE += -I$(MUTATIONPP_INSTALL)/include/mutation++
LIBS  = $(MUTATIONPP_INSTALL)/lib/libmutation++.so
EXEC  = \
	larsen

all: $(EXEC)

larsen:  larsen.cpp
	$(CC) -o larsen larsen.cpp $(INCLUDE) $(LIBS)
    
clean:
	rm -rf *.o
	rm -rf $(EXEC)

mrproper: clean
	rm -rf $(EXEC)
