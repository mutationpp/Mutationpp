
#include "mutation++.h"
#include "MultiPhaseEquilSolver.h"

using namespace Mutation::Thermodynamics;

int main()
{
    Mutation::Mixture mix("CO28");
    MultiPhaseEquilSolver equil(mix);
    
    double c[mix.nElements()];
    for (int i = 0; i < mix.nElements(); i++)
        c[i] = mix.getDefaultComposition(i);
    
    double x[mix.nSpecies()];
    std::pair<int, int> hist = equil.equilibrate(300.0, 101325.0, c, x);
    
    std::cout << "Steps: " << hist.first << std::endl;
    std::cout << "Newtons: " << hist.second << std::endl;
    
    std::cout << "Mole fractions: " << std::endl;
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << "X_" << mix.speciesName(i) << ": " << x[i] << std::endl;
}
