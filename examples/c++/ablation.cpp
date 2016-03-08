#include "mutation++.h"

int main() {

    Mutation::MixtureOptions opts("O2_ablation");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns   = mix.nSpecies();
    std::vector<double> v_rhoi(ns);
    std::vector<double> v_wall_production_rates(ns);
    std::vector<double> v_mole_fractions_edge(ns);

    double wall_temperature;

    const size_t state_var = 0;

    const int iO   = mix.speciesIndex("O");
    const int iO2  = mix.speciesIndex("O2");
    const int iCO  = mix.speciesIndex("CO");

// Conditions T = 3000K and p = 100Pa
    v_rhoi[iO]  = 6.00781e-05;
    v_rhoi[iO2] = 8.12943e-06;
    v_rhoi[iCO] = 1.00000e-07;
    double total_rho = v_rhoi[iO] + v_rhoi[iO2] + v_rhoi[iCO];
    wall_temperature = 1000.e0;
    

    // GET MOLE FRACTIONS
    int set_state_with_rhoi_T = 1;
    mix.setState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T);
    for (int i_ns ; i_ns < ns ; i_ns++){
        v_mole_fractions_edge[i_ns] = mix.X()[i_ns];
    }

    double gradient_distance_wall_bulk = 1.e0;
    mix.setWallState( &v_rhoi[0], &wall_temperature, state_var );
    mix.setDiffusionModel( &v_mole_fractions_edge[0], gradient_distance_wall_bulk );

    // SOLVE SURFACE BALANCE
    mix.solveSurfaceBalance();
    mix.getWallState( &v_rhoi[0], &wall_temperature, state_var );

    std::cout << "The result for the partial densities are: " << v_rhoi[0] << " " << v_rhoi[1] << std::endl;
    






















    mix.setWallState( &v_rhoi[0], &wall_temperature, state_var );
    mix.surfaceProductionRates( &v_wall_production_rates[0] );

    std::cout << v_wall_production_rates[0] << " " << v_wall_production_rates[1] << " " << v_wall_production_rates[2] << std::endl;

}
