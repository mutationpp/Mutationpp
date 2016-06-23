#include "mutation++.h"

int main() {

    Mutation::MixtureOptions opts( "N2_frc" );
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns   = mix.nSpecies();

    std::vector<double> v_rhoi(ns);
    std::vector<double> v_wall_prod_rates(ns);
    double wall_temperature;

    const int iN   = mix.speciesIndex("N");
    const int iN2  = mix.speciesIndex("N2");

// Conditions T = 3000K and p = 100Pa ??
    v_rhoi[iN]  = 6.00781e-05;
    v_rhoi[iN2] = 8.12943e-06;
    double total_rho = v_rhoi[iN] + v_rhoi[iN2];
    wall_temperature = 300.e0;
    
    int set_state_with_rhoi_T = 1;

    mix.setState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T); 
    mix.setWallState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T);

    mix.surfaceProductionRates(&v_wall_prod_rates[0]);
    std::cout << "The wall production rates are: " << v_wall_prod_rates[0] << " " << v_wall_prod_rates[1] << std::endl;

}
