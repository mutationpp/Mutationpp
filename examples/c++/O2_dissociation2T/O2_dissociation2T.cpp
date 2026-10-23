/**
 * @file O2_dissociation2T.cpp
 *
 * @brief O<sub>2</sub> Dissociation with 2T model example program.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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
    cout << setw(12) << "Tv[K]";
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
    cout << setw(12) << mix.Tv();
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
    const double T_init   = 10000.0;  // K
    const double Tv_init  = 300.0;  // K
    const double rho_init = 0.3899740209756568; // kg/m^3

    // First create the mixture object
    MixtureOptions opts;
    opts.setSpeciesDescriptor("O O2");        // 2 species O2 mixture
    opts.setThermodynamicDatabase("RRHO");    // Thermo database is RRHO
    opts.setMechanism("O2");                  // O2 + M = 2O + M
    opts.setStateModel("ChemNonEqTTv");       // chemical noneq. w/ 2T
    Mixture mix(opts);                        // Init. the mixture with opts

    // Setup arrays
    double rhoi [mix.nSpecies()];
    double wdot [mix.nSpecies()];
    double TToPass [2];
    double evib;
    double source;
    double EnergiesMass[2];

    // Set state of mixture to the initial conditions
    rhoi[mix.speciesIndex("O")]  = 0.0;
    rhoi[mix.speciesIndex("O2")] = 1.0*rho_init;

    TToPass[0] = T_init;
    TToPass[1] = Tv_init;
    mix.setState(rhoi, TToPass, 1); 
    mix.mixtureEnergies(EnergiesMass);

    EnergiesMass[0] = EnergiesMass[0]*mix.density();
    EnergiesMass[1] = EnergiesMass[1]*mix.density();

    // Write the results header and initial conditions
    printHeader(mix);
    printResults(0.0, mix);

    //time integration RK4
    double dt   = 1.0e-10;
    double tol  = 1.0e-14;
    double time = 0.0;
    double k1[mix.nSpecies()],k2[mix.nSpecies()],k3[mix.nSpecies()],k4[mix.nSpecies()], W[mix.nSpecies()], drhoi[mix.nSpecies()];
    double s1, s2, s3, s4, s[2], dEnergiesMass;
    double conv = 1.0;
    s[0] = EnergiesMass[0];
    s[1] = EnergiesMass[1];
    while (conv > tol) {

        dt = std::min(dt*1.0002, 1e-5);
        // Get the species production rates
        mix.setState(rhoi, EnergiesMass);
        mix.netProductionRates(k1);
        mix.energyTransferSource(&s1);

       for (int i = 0; i < mix.nSpecies(); ++i)
            W[i] = rhoi[i] + 0.5*dt*k1[i];
	s[1] = EnergiesMass[1] + 0.5*s1*dt;	

        mix.setState(W, s);
        mix.netProductionRates(k2);
        mix.energyTransferSource(&s2);

       for (int i = 0; i < mix.nSpecies(); ++i)
            W[i] = rhoi[i] + 0.5*dt*k2[i];
	s[1] = EnergiesMass[1] + 0.5*s2*dt;	

        mix.setState(W, s);
        mix.netProductionRates(k3);
        mix.energyTransferSource(&s3);

       for (int i = 0; i < mix.nSpecies(); ++i)
            W[i] = rhoi[i] + dt*k3[i];
	s[1] = EnergiesMass[1] + s3*dt;	

        mix.setState(W, s);
        mix.netProductionRates(k4);
        mix.energyTransferSource(&s4);

        conv = 0.0;
        for (int i = 0; i < mix.nSpecies(); ++i) {
            drhoi[i] = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0 * dt;
            rhoi[i] += drhoi[i];
            if (std::abs(drhoi[i]) > conv) conv = std::abs(drhoi[i]);
        }
        dEnergiesMass =  (s1 + 2.0*s2 + 2.0*s3 + s4) / 6.0 * dt;
	EnergiesMass[1] +=  dEnergiesMass;
        time += dt;

        // Set the new state of the mixture (using conserved variables)
        mix.setState(rhoi, EnergiesMass);
        printResults(time, mix);
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

