#include <iostream>

#include "mutation++.h"
using namespace Mutation;

void printModel(Mixture& mix)
{
    std::cout << "# of elements: " << mix.nElements() << '\n';
    std::cout << "# of species: " << mix.nSpecies() << ' ';
    std::cout << mix.nGas() << " (gas) " << mix.nCondensed() << " (condensed)\n";
    std::cout << "# of reactions: " << mix.nReactions() << '\n';
    std::cout << "# of temperatures: " << mix.nEnergyEqns() << "\n\n";

    std::cout << "Species:\n";
    for (auto& s: mix.species()) 
        std::cout << s.name() << ' ';
    std::cout << "\n\n";

    std::cout << "Reactions:\n";
    for (auto& r: mix.reactions())
        std::cout << r.formula() << '\n';
}

void printState(Mixture& mix)
{
    std::vector<double> T(mix.nEnergyEqns());
    mix.getTemperatures(&T[0]);

    std::cout << "\nTemperature(s) [K]: ";
    for (auto v: T) std::cout << v << ' ';
    
    std::cout << "\nPressure [Pa]: " << mix.P();

    std::cout << "\nSpecies mole fractions:\n";
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << mix.speciesName(i) << ": " << mix.X()[i] << '\n';
}

void printProperties(Mixture& mix)
{
    std::cout << "\nMixture enthalpy [J/kg]: " << mix.mixtureHMass();
    std::cout << "\nMixture entropy [J/kg-K]: " << mix.mixtureSMass();

    std::vector<double> lambda(mix.nEnergyEqns());
    mix.frozenThermalConductivityVector(&lambda[0]);
    std::cout << "\nMixture thermal conductivities [W/m-K]: ";
    for (auto v: lambda) std::cout << v << ' ';

    std::cout << "\nMixture dynamic viscosity [Pa-s]: " << mix.viscosity();
    std::cout << '\n';
}

void partOne()
{
    std::cout << "\n\n*** PART ONE ***\n";
    Mixture mix("air_5");
    printModel(mix);
    
    mix.equilibrate(300.0, ONEATM);
    printState(mix);
    printProperties(mix);

    mix.equilibrate(5000.0, ONEATM);
    printState(mix);
    printProperties(mix);

    std::vector<double> rho(mix.nSpecies()), T(mix.nEnergyEqns());
    mix.getTemperatures(&T[0]);
    mix.densities(&rho[0]);
    mix.setState(&rho[0], &T[0], 1);
    printState(mix);
    printProperties(mix);
}

void partTwo()
{
    std::cout << "\n\n*** PART TWO ***\n";
    MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setViscosityAlgorithm("Gupta-Yos");
    Mixture mix(opts);
    printModel(mix);

    mix.equilibrate(5000.0, ONEATM);
    printState(mix);
    printProperties(mix);
}

int main(int argc, char* argv[])
{
    partOne();
    partTwo();

    return 0;
}