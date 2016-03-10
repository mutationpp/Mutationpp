#include "mutation++.h"

int main() {


    std::cout << "@todo CHECK RATES OF NITRIDATION!" << std::endl;

    Mutation::MixtureOptions opts("O2_N2_ablation");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns   = mix.nSpecies();
    std::vector<double> v_rhoi(ns);
    std::vector<double> v_wall_production_rates(ns);
    std::vector<double> v_mole_fractions_edge(ns);

    double wall_temperature;

    std::vector<double> v_mole_fractions(ns);

    std::vector<double> v_diffusion_velocities(ns);
    std::vector<double> v_mole_gradients(ns);

    std::vector<double> v_function(ns);

    const size_t state_var = 0;

    const int iO   = mix.speciesIndex("O");
    const int iN  = mix.speciesIndex("N");
    const int iO2  = mix.speciesIndex("O2");
    const int iCO  = mix.speciesIndex("CO");
    const int iCN  = mix.speciesIndex("CN");

// Conditions T = 3000K and p = 100Pa
    v_rhoi[iO]  = 1.31228e-05;
    v_rhoi[iN] = 4.44957e-05;
    v_rhoi[iO2] = 3.87867e-07;
    v_rhoi[iCO] = 1.00000e-10;
    v_rhoi[iCN] = 1.00000e-10;
    double total_rho = v_rhoi[iO] + v_rhoi[iN] + v_rhoi[iO2] + v_rhoi[iCO] + v_rhoi[iCN];
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

    std::cout << "The result for the partial densities are: " << v_rhoi[0] << " " << v_rhoi[1] << " " << v_rhoi[2] 
              << " " << v_rhoi[3] << " " << v_rhoi[4] <<  std::endl;

//================================================================================================================

    total_rho = v_rhoi[iO] + v_rhoi[iN] + v_rhoi[iO2] + v_rhoi[iCO] + v_rhoi[iCN];

    mix.setState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T); 
    mix.setWallState(&v_rhoi[0], &wall_temperature, set_state_with_rhoi_T );

    for (int i_ns ; i_ns < ns ; i_ns++){
        v_mole_fractions[i_ns] = mix.X()[i_ns];
        v_mole_gradients[i_ns] = ( v_mole_fractions[i_ns] - v_mole_fractions_edge[i_ns] ) / gradient_distance_wall_bulk;
    }

    double electric_field = 0.E0;
    mix.stefanMaxwell(&v_mole_gradients[0], &v_diffusion_velocities[0], electric_field);
    mix.surfaceProductionRates( &v_wall_production_rates[0] );
    std::cout << "The wall production rates are: " << v_wall_production_rates[0] << " " 
              << v_wall_production_rates[1] << " " << v_wall_production_rates[2] << " " 
              << v_wall_production_rates[3] << " " << v_wall_production_rates[4] << std::endl;

    double mass_blowing = 0.0;
    for ( int i_ns = 0 ; i_ns < ns ; ++i_ns ){
        mass_blowing += v_wall_production_rates[i_ns];
    } 
    mass_blowing = mass_blowing / total_rho;

    for( int i_ns ; i_ns < ns ; i_ns++){
        v_function[i_ns] = v_rhoi[i_ns] * v_diffusion_velocities[i_ns] - v_wall_production_rates[i_ns]
                                                                       + mass_blowing * v_rhoi[i_ns]; 
    }

    std::cout << "The residual is equal to: " << v_function[0] << " " << v_function[1] 
              << " " << v_function[2] << " " << v_function[3] << " " << v_function[4] <<  std::endl;

}
