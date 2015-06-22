#include<cmath>
#include<mutation++.h>

/*
 * This file contains the implementation of some necessary functions for the test_gas_surface_interaction.cpp file
 */

/*
 * This function produces the pressure [10.0, 100000.0] and temperature [300, 3300] vectors for the test matrix
 */
void get_pressure_temperature_vectors ( std::vector<double>& v_pressure, std::vector<double>& v_temperature){
    const double minimum_pressure = log(10.0);
    const double maximum_pressure = log(100000.0);
    const double minimum_temperature = 300.0;
    const double maximum_temperature = 3000.0;

    for( size_t i_pressure = 0; i_pressure < v_pressure.size() ; ++i_pressure ){
        v_pressure[i_pressure] = std::exp( static_cast<double>(i_pressure) / static_cast<double>( v_pressure.size() - 1 ) * maximum_pressure + minimum_pressure );
    }
    for ( size_t i_temperature = 0; i_temperature < v_temperature.size(); ++i_temperature){
        v_temperature[i_temperature] = static_cast<double>(i_temperature) / static_cast<double>( v_temperature.size() - 1) * maximum_temperature + minimum_temperature;
    }

}

/*
 * This function computes the recombination rigorously for the the O-O2 recombination for the O2 mixture.
 */

double compute_production_rate_for_O_O2_recombination( const double& gamma, const std::vector<double>& Mw_O, const std::vector<double>& rho_O, const double& T ){
    int position_of_oxygen_in_mixture_file = 0;
    return gamma * rho_O[ position_of_oxygen_in_mixture_file ] * sqrt( Mutation::RU / Mw_O[ position_of_oxygen_in_mixture_file ]  * T  / ( 2.0 * Mutation::PI ) );
}

/*
 * This function computes the recombination rigorously for the the N-N2 recombination for the air5 mixture.
 */
double compute_production_rate_for_N_N2_recombination_air( const double& gamma_NN, const std::vector<double>& Mw_N, const std::vector<double>& rho_N, const double& T ){
    int position_of_nitrogen_in_mixture_file = 0;

    return gamma_NN * rho_N[ position_of_nitrogen_in_mixture_file ] * sqrt( Mutation::RU / Mw_N[ position_of_nitrogen_in_mixture_file ]  * T  / ( 2.E0 * Mutation::PI ) );
}

/*
 * This function computes the recombination rigorously for the the O-O2 recombination for the air5 mixture.
 */
double compute_production_rate_for_O_O2_recombination_air( const double& gamma_OO, const std::vector<double>& Mw_O, const std::vector<double>& rho_O, const double& T ){
    int position_of_oxygen_in_mixture_file = 1;

    return gamma_OO * rho_O[ position_of_oxygen_in_mixture_file ] * sqrt( Mutation::RU / Mw_O[ position_of_oxygen_in_mixture_file ]  * T  / ( 2.E0 * Mutation::PI ) );
}

/*
 * This function sums over all species the production rates.
 */
double sum_production_rates( const std::vector<double>& p_production_rates ){
    double sum = 0.E0;
    for (int i_sp = 0; i_sp < p_production_rates.size(); i_sp++ ) {
        sum += p_production_rates[i_sp];
    }
    return sum;
}
