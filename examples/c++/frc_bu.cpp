
#include "mutation++.h"
#include <iostream>
#include <cmath>

int main(){

    Mutation::MixtureOptions opts("N2_frc");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();

    std::vector<double> v_rhoi(ns);
    double wall_temperature;

    const int iN = mix.speciesIndex("N");
    const int iN2 = mix.speciesIndex("N2");

    v_rhoi[iN]  = 6.00781e-05;
    v_rhoi[iN2] = 8.12943e-06;
    double total_rho = v_rhoi[iN] + v_rhoi[iN2];
    wall_temperature = 300.e0;

    int set_state_with_rhoi_T = 1;

//================================================================================================================


    return 0;

}


