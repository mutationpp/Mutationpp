#include "mutation++.h"

#include <iostream>

int main()
{
   /* Mutation::MixtureOptions opts( "Meteor_wall" );
    opts.setStateModel("EquilTP");
    opts.setThermodynamicDatabase("NASA-9");

    Mutation::Mixture* mix = new Mutation::Mixture(opts);*/

    Mutation::Mixture mix("Carbon");
    const int ns   = mix.nSpecies();


    double wall_temperature = 1500.E0;
    double pressure = 0.01*Mutation::ONEATM;
    int state = 0;

    double BprimeChar;
    std::vector<double> v_rhoi(ns);
    std::vector<double> v_xi(ns);
    std::vector<double> SpeciesCharMass;
    std::vector<std::string> SpeciesCharName;

    for(int i = 0; i < ns; ++i){

        v_rhoi[i] = 0.E0;
    }



    mix.setWallState(&pressure, &wall_temperature, state);
    mix.solveSurfaceBalance();
    mix.getWallState(&v_rhoi[0], &wall_temperature, 1);
    mix.BprimeCharSpecies(SpeciesCharName);
    mix.BprimeSolution(BprimeChar, SpeciesCharMass);
    //mix.BprimeSolution(BprimeChar,SpeciesCharName, SpeciesCharMass);

    mix.Thermodynamics::convert<Mutation::Thermodynamics::RHO_TO_X>(&v_rhoi[0],&v_xi[0]);

    std::cout << std::setw(10) << "\"Tw[K]\""
         << std::setw(15) << "\"B'c\"";
    for (int i = 0; i < ns; ++i){
        std::cout << std::setw(25) << "\"" + mix.speciesName(i) + "\"";}
    for (int i = 0; i < SpeciesCharName.size(); ++i){
        std::cout << std::setw(25) << "\"" + SpeciesCharName[i] + "\"";}
    std::cout << std::endl;

    std::cout << std::setw(10) << wall_temperature << std::setw(15) << BprimeChar;
    for (int i = 0; i < ns; ++i){
        std::cout << std::setw(25) << v_xi[i];}
    for (int i = 0; i < SpeciesCharMass.size(); ++i){
        std::cout << std::setw(25) << SpeciesCharMass[i];}
    std::cout << std::endl;
}



