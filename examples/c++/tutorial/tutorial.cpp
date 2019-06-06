#include <iostream>

#include "mutation++.h"
using namespace Mutation;

int main(int argc, char* argv[])
{
    Mixture mix("air_5");

    std::cout << "# of elements: " << mix.nElements() << '\n';
    std::cout << "# of species: " << mix.nSpecies() << ' ';
    std::cout << mix.nGas() << " (gas) " << mix.nCondensed() << " (condensed)\n";
    std::cout << "# of reactions: " << mix.nReactions() << '\n';
    std::cout << "# of temperatures: " << mix.nEnergyEqns() << '\n';

    std::cout << "Species:\n";
    for (auto& s: mix.species()) 
        std::cout << s.name() << ' ';
    std::cout << "\n\n";

    std::cout << "Reactions:\n";
    for (auto& r: mix.reactions())
        std::cout << r.formula() << '\n';
    std::cout << '\n';
    
    
    return 0;
}