#include "mutation++.h"

int main() {

    Mutation::MixtureOptions opts( "O2_gamma_test" );
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns   = mix.nSpecies();
    std::vector<double> v_rhoi(ns);
    std::vector<double> v_wall_production_rates(ns);
    std::vector<double> v_mole_fractions(ns);
    std::vector<double> v_mole_fractions_edge(ns);

    std::vector<double> v_diffusion_velocities(ns);
    std::vector<double> v_mole_gradients(ns);

    std::vector<double> v_function(ns);

    double bulk_temperature;
    double wall_temperature;

    const int iO   = mix.speciesIndex("O");
    const int iO2  = mix.speciesIndex("O2");

// Conditions T = 3000K and p = 100Pa
    v_rhoi[iO]  = 6.00781e-05;
    v_rhoi[iO2] = 8.12943e-06;
    double total_rho = v_rhoi[iO] + v_rhoi[iO2];
    wall_temperature = 300.e0;
    
    int set_state_with_rhoi_T = 1;
    mix.setState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T);
    for (int i_ns ; i_ns < ns ; i_ns++){
        v_mole_fractions_edge[i_ns] = mix.X()[i_ns];
    }

    int state_var = 0;
    double gradient_distance_wall_bulk = 1.e0;
    mix.setWallState( &v_rhoi[0], &wall_temperature, state_var );
    mix.setDiffusionModel( &v_mole_fractions_edge[0], gradient_distance_wall_bulk );

//    std::cout << v_wall_production_rates[0] << " " << v_wall_production_rates[1] << " CHECK ME" << std::endl;


}
