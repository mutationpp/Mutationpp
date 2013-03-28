
#include "mutation++.h"
#include "MultiPhaseEquilSolver.h"

using namespace Mutation::Thermodynamics;

int main()
{
    std::vector<std::string> species;
    species.push_back("CO");
    species.push_back("CO2");
    species.push_back("O2");
    species.push_back("C(s)");
    
    Thermodynamics thermo(species, "NASA-9", "T");
    MultiPhaseEquilSolver equil(thermo);
    
    double c[2];
    double x[4];
    equil.equilibrate(3000.0, 101325.0, c, x);
}
