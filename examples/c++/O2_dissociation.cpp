/**
 * \example O2_dissociation.cpp
 *
 * Example program demonstrating the use of the Mutation++ finite-rate 
 * chemistry routine netProductionRates().  This program simulates
 * a constant volume, adiabatic reactor beginning with pure O2 at 4000 K and
 * density of 2.85e-4 kg/m^3.  The program integrates the mass conservation
 * equation for each species explicitly in time up to 100 seconds until the O/O2
 * mixture is essentially at equilibrium.  The reactor temperature, pressure,
 * and mass fractions of O and O2 are output to show the results of the
 * simulation.  The final conditions at 100 seconds are also compared with
 * equilibrium values for the O/O2 mixture at the final temperature and
 * pressure.
 *
 * The main example code is illustrated as follows:
 *
 * Running the example yields the following output:
 */

// Must include this header file to use the Mutation++ library
#include "mutation++.h"

using namespace Mutation;
using namespace std;

// Prints the header of the output table
void printHeader(const Mixture& mix)
{
    cout << setw(15) << "Time[s]";
    cout << setw(15) << "T[K]";
    cout << setw(15) << "P[Pa]";
    cout << setw(15) << "rho[kg/m^3]";
    cout << setw(15) << "e[J/kg]";
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(15) << "Y_" + mix.speciesName(i);
    cout << endl;
}

// Prints the mixture properties at the given time
void printResults(double time, const Mixture& mix)
{
    cout << setw(15) << time;
    cout << setw(15) << mix.T();
    cout << setw(15) << mix.P();
    cout << setw(15) << mix.density();
    cout << setw(15) << mix.mixtureEnergyMass();
    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(15) << mix.Y()[i];
    cout << endl;
}

//! [main]
int main()
{
    // Initial conditions are defined here
    const double T_init   = 4000.0;  // K
    const double rho_init = 2.85e-4; // kg/m^3

    // First create the mixture object
    MixtureOptions opts;
    opts.setSpeciesDescriptor("O O2");        // 2 species O2 mixture
    opts.setThermodynamicDatabase("NASA-9");  // NASA-9 coeff. poly. thermo
    opts.setMechanism("O2");                  // O2 + M = 2O + M
    opts.setStateModel("ChemNonEq1T");        // chemical noneq. w/ 1T
    Mixture mix(opts);                        // Init. the mixture with opts

    // Setup arrays
    double rhoi [mix.nSpecies()];
    double wdot [mix.nSpecies()];

    // Set state of mixture to the initial conditions
    rhoi[mix.speciesIndex("O")]  = 0.0;
    rhoi[mix.speciesIndex("O2")] = rho_init;
    mix.setState(rhoi, &T_init, 1); // set state using {rho_i, T}

    // Get the mixture energy which is conserved during the simulation
    double etot = mix.mixtureEnergyMass() * mix.density(); // J/m^3

    // Write the results header and initial conditions
    printHeader(mix);
    printResults(0.0, mix);

    // Integrate the species mass conservation equations forward in time
    double dt   = 0.0;
    double time = 0.0;
    double tout = 0.0001;

    while (time < 100.0) {
        // Get the species production rates
        mix.netProductionRates(wdot);

        // Compute "stable" time-step based on maximum allowed change in
        // species densities
        dt = 0.0;
        for (int i = 0; i < mix.nSpecies(); ++i)
            dt += 1.0e-7 * std::abs(rho_init * mix.Y()[i] / wdot[i]);
        dt = min(dt / mix.nSpecies(), tout-time);

        // Integrate in time
        for (int i = 0; i < mix.nSpecies(); ++i)
            rhoi[i] += wdot[i]*dt;
        time += dt;

        // Set the new state of the mixture (using conserved variables)
        mix.setState(rhoi, &etot);

        // Print results at specified time intervals
        if (std::abs(tout-time) < 1.0e-10) {
           printResults(time, mix);

           if (tout > 9.9)          tout += 10.0;
           else if (tout > 0.99)    tout += 1.0;
           else if (tout > 0.099)   tout += 0.1;
           else if (tout > 0.0099)  tout += 0.01;
           else if (tout > 0.00099) tout += 0.001;
           else                     tout += 0.0001;
        }
    }

    // Now get the equilibrium values
    mix.equilibriumComposition(mix.T(), mix.P(), rhoi); // puts X_i in rhoi
    mix.convert<Thermodynamics::X_TO_Y>(rhoi, rhoi);    // converts X_i to Y_i

    cout << "Equilibrium mass fractions at " << mix.T() << " K and "
         << mix.P() << " Pa:" << endl;

    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(15) << mix.speciesName(i) << ":"
             << setw(15) << rhoi[i] << endl;

    return 0;
}
//! [main]

