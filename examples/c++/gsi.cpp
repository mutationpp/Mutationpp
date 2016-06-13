#include "mutation++.h"

int main() {

    Mutation::MixtureOptions opts( "O2_gamma_test" );
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns   = mix.nSpecies();
    Eigen::VectorXd v_rhoi(ns);
    Eigen::VectorXd v_wall_production_rates(ns);
    Eigen::VectorXd v_mole_fractions(ns);
    Eigen::VectorXd v_mole_fractions_edge(ns);

    Eigen::VectorXd v_diffusion_velocities(ns);
    Eigen::VectorXd v_mole_gradients(ns);

    Eigen::VectorXd v_function(ns);

    double bulk_temperature;
    double wall_temperature;

    const int iO   = mix.speciesIndex("O");
    const int iO2  = mix.speciesIndex("O2");

// Conditions T = 3000K and p = 100Pa
    v_rhoi(iO)  = 6.00781e-05;
    v_rhoi(iO2) = 8.12943e-06;
    double total_rho = v_rhoi.sum();
    wall_temperature = 300.e0;
    
    int set_state_with_rhoi_T = 1;
    mix.setState( v_rhoi.data(), &wall_temperature, set_state_with_rhoi_T);
    for (int i_ns = 0 ; i_ns < ns ; i_ns++){
        v_mole_fractions_edge(i_ns) = mix.X()[i_ns];
    }

    int state_var = 0;
    double gradient_distance_wall_bulk = 1.e0;
    mix.setWallState( v_rhoi.data(), &wall_temperature, state_var );
    mix.setDiffusionModel( v_mole_fractions_edge.data(), gradient_distance_wall_bulk );

    mix.solveSurfaceBalance();
    mix.getWallState( v_rhoi.data(), &wall_temperature, state_var );

    std::cout << "The result for the partial densities are: " << v_rhoi(0) << " " << v_rhoi(1) << std::endl;

//================================================================================================================

    mix.setState( v_rhoi.data(), &wall_temperature, set_state_with_rhoi_T); 
    mix.setWallState( v_rhoi.data(), &wall_temperature, set_state_with_rhoi_T );

    v_mole_fractions = Eigen::Map<const Eigen::VectorXd>(mix.X(), ns);
    v_mole_gradients = ( v_mole_fractions - v_mole_fractions_edge ) / gradient_distance_wall_bulk;

    double electric_field = 0.E0;
    mix.stefanMaxwell( v_mole_gradients.data(), v_diffusion_velocities.data(), electric_field);
    v_wall_production_rates = Eigen::Map<const Eigen::VectorXd>(mix.surfaceProductionRates(), ns);
    std::cout << "The wall production rates are: " << v_wall_production_rates(0)
                                            << " " << v_wall_production_rates(1) << std::endl;

    v_function = v_rhoi.cwiseProduct( v_diffusion_velocities ) - v_wall_production_rates;

    std::cout << "The residual is equal to: "<< v_function(0) << " " << v_function(1) << std::endl;

//================================================================================================================

}
