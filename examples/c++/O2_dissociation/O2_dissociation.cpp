/**
 * @file O2_dissociation.cpp
 *
 * @brief O<sub>2</sub> Dissociation example program.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/**
 * @page example_O2_dissociation O<sub>2</sub> Dissociation
 *
 * @tableofcontents
 *
 * @section O2_diss_example_intro Introduction
 * This example program simulates a constant volume, adiabatic reactor beginning
 * with pure O<sub>2</sub> at 4000 K and density of 2.85e-4 kg/m^3.  The program
 * integrates the mass conservation equation for each species explicitly in time
 * up to 100 seconds until the O/O<sub>2</sub> mixture is essentially at
 * equilibrium.  The reactor temperature, pressure, and mass fractions of O and
 * O<sub>2</sub> are output to show the results of the simulation.  The final
 * conditions at 100 seconds are also compared with equilibrium values for the
 * O/O<sub>2</sub> mixture at the final temperature and pressure.
 *
 * The example demonstrates the following features:
 * - Using the Mutation::MixtureOptions class to modify the basic options in a
 *   [mixture](@ref mixtures)
 * - Using the [ChemNonEq1T](@ref Mutation::Thermodynamics::ChemNonEqStateModel)
 *   state model for chemical nonequilibrium
 * - Using Mutation::Mixture::netProductionRates() to compute mass production
 *   rates of each species due to chemical reactions
 * - Setting the state of a mixture
 * - Retrieving various properties of a mixture
 * - Using Mutation::Mixture::convert() to convert species mole to mass
 *   fractions
 *
 * @section O2_diss_example_code Example Code
 * @snippet examples/c++/O2_dissociation.cpp example_code
 *
@section O2_diss_example_output Expected Output
@verbatim
 Time[s]        T[K]       P[Pa] rho[kg/m^3]     e[J/kg]         Y_O        Y_O2
       0        4000     296.214    0.000285 3.22221e+06           0           1
  0.0001     3969.25     294.531    0.000285 3.22221e+06  0.00201901    0.997981
  0.0002     3941.29     292.991    0.000285 3.22221e+06  0.00385366    0.996146
  0.0003     3915.68     291.574    0.000285 3.22221e+06  0.00553232    0.994468
  0.0004      3892.1     290.263    0.000285 3.22221e+06  0.00707747    0.992923
  0.0005     3870.26     289.045    0.000285 3.22221e+06  0.00850721    0.991493
  0.0006     3849.95     287.907    0.000285 3.22221e+06   0.0098363    0.990164
  0.0007     3830.98      286.84    0.000285 3.22221e+06   0.0110769    0.988923
  0.0008     3813.14     285.833    0.000285 3.22221e+06   0.0122392    0.987761
  0.0009     3796.34      284.88    0.000285 3.22221e+06   0.0133314    0.986669
   0.001     3780.49     283.979    0.000285 3.22221e+06   0.0143608    0.985639
   0.002      3658.8     276.975    0.000285 3.22221e+06   0.0222508    0.977749
   0.003      3576.3     272.138    0.000285 3.22221e+06    0.027563    0.972437
   0.004      3514.8     268.485    0.000285 3.22221e+06    0.031513    0.968487
   0.005     3466.01     265.559    0.000285 3.22221e+06   0.0346303     0.96537
   0.006      3425.9     263.135    0.000285 3.22221e+06   0.0371893    0.962811
   0.007     3391.95     261.071    0.000285 3.22221e+06   0.0393508    0.960649
   0.008     3362.57     259.273    0.000285 3.22221e+06   0.0412152    0.958785
   0.009     3336.77     257.688    0.000285 3.22221e+06   0.0428501     0.95715
    0.01     3313.83     256.273    0.000285 3.22221e+06   0.0443031    0.955697
    0.02     3167.01     247.085    0.000285 3.22221e+06   0.0535421    0.946458
    0.03     3085.68     241.902    0.000285 3.22221e+06   0.0586269    0.941373
    0.04     3030.15     238.323    0.000285 3.22221e+06   0.0620773    0.937923
    0.05     2988.47     235.617    0.000285 3.22221e+06   0.0646619    0.935338
    0.06     2955.26     233.448    0.000285 3.22221e+06   0.0667151    0.933285
    0.07     2927.77     231.644    0.000285 3.22221e+06   0.0684097     0.93159
    0.08     2904.43     230.106    0.000285 3.22221e+06   0.0698477    0.930152
    0.09     2884.19     228.769    0.000285 3.22221e+06   0.0710933    0.928907
     0.1     2866.32     227.584    0.000285 3.22221e+06   0.0721896     0.92781
     0.2     2755.06     220.137    0.000285 3.22221e+06   0.0789885    0.921011
     0.3     2694.92      216.06    0.000285 3.22221e+06     0.08264     0.91736
     0.4      2654.5     213.301    0.000285 3.22221e+06   0.0850858    0.914914
     0.5     2624.49     211.241    0.000285 3.22221e+06    0.086897    0.913103
     0.6      2600.9     209.616    0.000285 3.22221e+06   0.0883192    0.911681
     0.7     2581.59     208.282    0.000285 3.22221e+06   0.0894791    0.910521
     0.8     2565.41     207.161    0.000285 3.22221e+06   0.0904507    0.909549
     0.9     2551.57     206.201    0.000285 3.22221e+06    0.091281    0.908719
       1     2539.56     205.366    0.000285 3.22221e+06   0.0920014    0.907999
       2     2470.36     200.526    0.000285 3.22221e+06   0.0961351    0.903865
       3     2439.67     198.365    0.000285 3.22221e+06   0.0979612    0.902039
       4     2423.24     197.204    0.000285 3.22221e+06   0.0989377    0.901062
       5     2413.73     196.531    0.000285 3.22221e+06   0.0995024    0.900498
       6     2408.01     196.126    0.000285 3.22221e+06   0.0998421    0.900158
       7     2404.49     195.876    0.000285 3.22221e+06    0.100051    0.899949
       8     2402.29      195.72    0.000285 3.22221e+06    0.100181    0.899819
       9     2400.91     195.623    0.000285 3.22221e+06    0.100263    0.899737
      10     2400.04     195.561    0.000285 3.22221e+06    0.100315    0.899685
      20     2398.55     195.455    0.000285 3.22221e+06    0.100403    0.899597
      30     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
      40     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
      50     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
      60     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
      70     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
      80     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
      90     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
     100     2398.54     195.454    0.000285 3.22221e+06    0.100404    0.899596
Equilibrium mass fractions at 2398.54 K and 195.454 Pa:
              O:       0.100403
             O2:       0.899597
@endverbatim
 */

/// [example_code]
// Must include this header file to use the Mutation++ library
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
    MixtureOptions opts;
    opts.setSpeciesDescriptor("O O2");        // 2 species O2 mixture
    opts.setThermodynamicDatabase("RRHO");    // Thermo database is RRHO
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
    mix.convert<X_TO_Y>(rhoi, rhoi);                    // converts X_i to Y_i

    cout << "Equilibrium mass fractions at " << mix.T() << " K and "
         << mix.P() << " Pa:" << endl;

    for (int i = 0; i < mix.nSpecies(); ++i)
        cout << setw(15) << mix.speciesName(i) << ":"
             << setw(15) << rhoi[i] << endl;

    return 0;
}
/// [example_code]

