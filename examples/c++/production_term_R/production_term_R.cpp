#include "mutation++.h"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace std;

// Prints the header of the output table
void printHeader(const Mixture& mix)
{
    cout << setw(8) << "Time[s]";
    cout << setw(12) << "T[K]";
    cout << setw(12) << "P[Pa]";
    cout << setw(12) << "rho[kg/m^3]";
    cout << setw(12) << "e[J/kg]";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << "Y_" + mix.speciesName(i);
    cout << endl;
}

// Prints the mixture properties at the given time
void printResults(double time, const Mixture& mix)
{
    cout << setw(8) << time;
    cout << setw(12) << mix.T();
    cout << setw(12) << mix.P();
    cout << setw(12) << mix.density();
    cout << setw(12) << mix.mixtureEnergyMass();
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(12) << mix.Y()[i];
    cout << endl;
}

// Main entry point
int main()
{
    // Initial conditions are defined here
    const double T_init   = 4000.0;  // K
    const double rho_init = 2.85e-4; // kg/m^3

    // First create the mixture object
    MixtureOptions opts("air_5");
    opts.setThermodynamicDatabase("RRHO_3Tv");
    opts.setMechanism("air5");                 
    opts.setStateModel("ChemNonEqT3Tv");        
    Mixture mix(opts);                       

    // Setup arrays
    double rhoi [mix.nSpecies()];
    //double wdot [mix.nSpecies()];
    double wdot [mix.nMolecules()];

    double Ts[4];
    Ts[0] = 5000; 
    Ts[1] = 6000;
    Ts[2] = 6000; 
    Ts[3] = 6000; 

    // Set state of mixture to the initial conditions
    rhoi[mix.speciesIndex("N")]  = 0.20;
    rhoi[mix.speciesIndex("O")]  = 0.20;
    rhoi[mix.speciesIndex("NO")]  = 0.20;
    rhoi[mix.speciesIndex("N2")]  = 0.20;
    rhoi[mix.speciesIndex("O2")]  = 0.20;

    mix.setState(rhoi, Ts, 1); 

    // Get the mixture energy which is conserved during the simulation
    //double etot = mix.mixtureEnergyMass() * mix.density(); // J/m^3

    // Write the results header and initial conditions
    printHeader(mix);
    printResults(0.0, mix);

    mix.R(wdot);
    //mix.R_22(wdot);
    //mix.R_23(wdot);

    for (int i=0; i<mix.nMolecules(); ++i)
        std::cout <<  " R["<<i<<"] = " << wdot[i] << std::endl;

    return 0;
}
