
#include "mutation++.h"

#include<memory>
#include<iostream>
#include<vector>
#include<cstdio>
#include<cassert>
#include<fstream>
#include<functional>

using namespace std;
using namespace Mutation;
using namespace Mutation::Thermodynamics;

// Main entry point
int main()
{
    // First create the mixture object
    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2(*)");		
    opts.setThermodynamicDatabase("RRaHO");    
    // opts.setMechanism("N2");               
    // opts.setStateModel("ChemNonEqvSTS");    
    Mixture mix(opts);                     

    double P = ONEATM;
    double T = 300.0;
    std::cout << std::setw(7) << mix.T() << std::endl;
    mix.setState(&T, &P);

    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(15) << mix.speciesName(i) << " " 
             << setw(15) << mix.Y()[i] << " " 
             << setw(15) << mix.speciesMw(i) << " " 
	     << endl;

    cout << mix.nSpecies() << endl;

    return 0;
}
