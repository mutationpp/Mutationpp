
#include "mutation++.h"
#include "MultiPhaseEquilSolver.h"

using namespace Mutation::Thermodynamics;
using namespace std;

int main()
{
    Mutation::Mixture mix("CO28");
    MultiPhaseEquilSolver equil(mix);
    
    double c[mix.nElements()];
    for (int i = 0; i < mix.nElements(); i++)
        c[i] = mix.getDefaultComposition(i);
    
    double x[mix.nSpecies()];
    std::pair<int, int> hist = equil.equilibrate(5000.0, 101325.0, c, x);
    
    std::cout << "Steps: " << hist.first << std::endl;
    std::cout << "Newtons: " << hist.second << std::endl;
    
    std::cout << "Mole fractions: " << std::endl;
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << "X_" << mix.speciesName(i) << ": " << x[i] << std::endl;
        
    std::cout << std::endl;
    std::cout << "With constraint X_C = X_O" << std::endl;
    
    double A[mix.nSpecies()];
    for (int i = 0; i < mix.nSpecies(); ++i)
        A[i] = 0.0;
    A[mix.speciesIndex("CO2")] = 1.0;
    A[mix.speciesIndex("O")] = -1.0;
    
    equil.addConstraint(A, 0.0);
    
    hist = equil.equilibrate(5000.0, 101325.0, c, x);
    
    std::cout << "Steps: " << hist.first << std::endl;
    std::cout << "Newtons: " << hist.second << std::endl;
    
    std::cout << "Mole fractions: " << std::endl;
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << "X_" << mix.speciesName(i) << ": " << x[i] << std::endl;
}
