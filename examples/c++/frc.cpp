#include "mutation++.h"

int main() {

    Mutation::MixtureOptions opts( "N2_frc" );
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);


    const int ns   = mix.nSpecies();


    std::vector<double> v_rhoi(ns);
/*
    std::vector<double> v_wall_production_rates(ns);
    std::vector<double> v_mole_fractions(ns);
    std::vector<double> v_mole_fractions_edge(ns);

    std::vector<double> v_diffusion_velocities(ns);
    std::vector<double> v_mole_gradients(ns);

    std::vector<double> v_function(ns);

    double bulk_temperature;
*/
    double wall_temperature;


    const int iN   = mix.speciesIndex("N");
    const int iN2  = mix.speciesIndex("N2");

// Conditions T = 3000K and p = 100Pa ??
    v_rhoi[iN]  = 0.0e0; //6.00781e-05;
    v_rhoi[iN2] = 0.0e0;// 8.12943e-06;
    double total_rho = v_rhoi[iN] + v_rhoi[iN2];
    wall_temperature = 300.e0;
    
    int set_state_with_rhoi_T = 1;
//    mix.setState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T);
//    for (int i_ns ; i_ns < ns ; i_ns++){
//        v_mole_fractions_edge[i_ns] = mix.X()[i_ns];
//    }

/*
    int state_var = 0;
    double gradient_distance_wall_bulk = 1.e0;
    mix.setWallState( &v_rhoi[0], &wall_temperature, state_var );
    mix.setDiffusionModel( &v_mole_fractions_edge[0], gradient_distance_wall_bulk );

    mix.solveSurfaceBalance();
    mix.getWallState( &v_rhoi[0], &wall_temperature, state_var );

    std::cout << "The result for the partial densities are: " << v_rhoi[0] << " " << v_rhoi[1] << std::endl;

//================================================================================================================

    mix.setState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T); 
    mix.setWallState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T );

    for (int i_ns ; i_ns < ns ; i_ns++){
        v_mole_fractions[i_ns] = mix.X()[i_ns];
        v_mole_gradients[i_ns] = ( v_mole_fractions[i_ns] - v_mole_fractions_edge[i_ns] ) / gradient_distance_wall_bulk;
    }

    double electric_field = 0.E0;
    mix.stefanMaxwell(&v_mole_gradients[0], &v_diffusion_velocities[0], electric_field);
    mix.surfaceProductionRates( &v_wall_production_rates[0] );
    std::cout << "The wall production rates are: " << v_wall_production_rates[0] << " " << v_wall_production_rates[1] << std::endl;

    for( int i_ns ; i_ns < ns ; i_ns++){
        v_function[i_ns] = v_rhoi[i_ns] * v_diffusion_velocities[i_ns] - v_wall_production_rates[i_ns];
    }

    std::cout << "The residual is equal to: "<< v_function[0] << " " << v_function[1] << std::endl;
*/

}
