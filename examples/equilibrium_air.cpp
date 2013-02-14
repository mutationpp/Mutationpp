/**
 * @file equilibrium_air.cpp
 *
 * Example program which makes use of the Mutation++ library to compute
 * equilibrium mole fractions, mixture frozen specific heat, enthalpy, and 
 * entropy for an air mixture containing 20% O and 80% N at a pressure of 1 atm
 * and a range of temperatures from 300 to 15,000 K using the NASA 9-coefficient 
 * polynomial thermodynamic database.  The example demonstrates using the 
 * MixtureOptions feature to change the default library options and how to 
 * compute properties for an equilibrium mixture.  The output is a property
 * table written to standard output as a function of temperature.
 */


// Must include this header file to use the Mutation++ library
#include "mutation++.h"
#include <iostream>

int main()
{
    // Generate the default options for the air11 mixture
    MixtureOptions opts("air11");
    
    // Override composition in mixture file
    opts.setDefaultComposition()
        ("e-", 0.0)
        ("N",  0.8)
        ("O",  0.2);
    
    // Change from default thermodynamic database (RRHO) to NASA 9-coefficient
    // polynomial database
    opts.setThermodynamicDatabase("NASA-9");
    
    // Load the mixture with the new options
    Mixture mix(opts);
    
    // Write a header line for the table
    std::cout << std::setw(15) << "T[K]";
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << std::setw(15) << "X_" + mix.speciesName(i);
    std::cout << std::setw(15) << "Cp,fr[J/kg-K]";
    std::cout << std::setw(15) << "H[J/kg]";
    std::cout << std::setw(15) << "S[J/kg-K]";
    std::cout << std::endl;
    
    // Loop over range of temperatures and compute equilibrium values at 1 atm    
    for (int i = 0; i < 148; ++i) {
        // Compute the temperature
        double T = 300.0 + static_cast<double>(i) * 100.0;
        
        // Set the mixture state equal to the equilibrium state for the given
        // temprature and pressure
        mix.equilibrate(T, ONEATM);
        
        // Temperature
        std::cout << std::setw(15) << mix.T();
        
        // Species mole fractions
        for (int j = 0; j < mix.nSpecies(); ++j)
            cout << std::setw(15) << mix.X()[j];
        
        // Other properties
        std::cout << std::setw(15) << mix.mixtureFrozenCpMass(); // Cp [J/kg-K]
        std::cout << std::setw(15) << mix.mixtureHMass();        // H  [J/kg]
        std::cout << std::setw(15) << mix.mixtureSMass();        // S  [J/kg-K]
        
        std::cout << std::endl;
    }
    
    return 0;
}

