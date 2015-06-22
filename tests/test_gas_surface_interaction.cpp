
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE catalysis tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <vector>
// #include <numeric>
#include "MixtureFixtures.h"
#include "CatalysisFunctions.h"

BOOST_AUTO_TEST_SUITE(gsi_catalysis)

//  BOOST_FIXTURE_TEST_CASE(gsi_reactions_O2_gamma, SetMixtureGSIO2gamma)
//   {
//       const double tol = 1.0e3*std::numeric_limits<double>::epsilon();
//   
//       const int nb_ns = mix().nSpecies();
//       const int nt = mix().nEnergyEqns();
//       std::vector<double> rhoi(nb_ns);
//       std::vector<double> Mw(nb_ns);
//       std::vector<double> p_production_rates(nb_ns);
//       double sum = 0.0;
//       
//       for (int i_sp = 0; i_sp < nb_ns; i_sp++){ 
//           Mw[i_sp] = mix().speciesMw(i_sp);
//       }
//   
//       int nb_gsi_reactions = mix().nCatalyticReactions();
//       BOOST_CHECK_EQUAL(nb_gsi_reactions, 1); 
//       double gamma_coefficient = 0.01;
//       BOOST_CHECK_CLOSE(gamma_coefficient, 0.01, tol);
//   
//       int nb_pressures = 10;
//       int nb_temperatures = 10;
//       std::vector<double> v_pressure( nb_pressures );
//       std::vector<double> v_temperature( nb_temperatures );
//       get_pressure_temperature_vectors( v_pressure, v_temperature );
//   
//       for (int i_pressure = 0 ; i_pressure < nb_pressures; ++i_pressure) {
//           for (int i_temperature = 0 ; i_temperature < nb_temperatures; ++i_temperature) {
//   
//               equilibrateMixture(v_temperature[i_temperature], v_pressure[i_pressure]);
//               rhoi = getPartialDensities();
//               mix().setWallState(&rhoi[0],&v_temperature[i_temperature]);
//               mix().netGSIProductionRates( &p_production_rates[0] );
//               sum = sum_production_rates( p_production_rates );
//   
//   //            std::cout << "The pressure is equal to p = " << v_pressure[i_pressure] << " and the temperature equal to T = " << v_temperature[i_temperature] << std::endl;
//   //            std::cout << "The rhois are equal to rhoi[0] = " << rhoi[0] << " and rhoi[1] = " << rhoi[1] << std::endl;
//   
//               BOOST_CHECK_CLOSE(p_production_rates[0], compute_production_rate_for_O_O2_recombination(gamma_coefficient, Mw, rhoi, v_temperature[i_temperature]), tol );
//               BOOST_CHECK_CLOSE(sum, 0.0, tol);
//   
//           }
//       }
//   
//   }

// BOOST_FIXTURE_TEST_CASE(gsi_reactions_air5_gamma, SetMixtureGSIair5gamma)
// {
//     const double tol = 1.0e3*std::numeric_limits<double>::epsilon();
// 
//     const int nb_ns = mix().nSpecies();
//     const int nt = mix().nEnergyEqns();
//     std::vector<double> rhoi(nb_ns);
//     std::vector<double> Mw(nb_ns);
//     std::vector<double> p_production_rates(nb_ns);
//     double sum;
//     
//     for (int i_sp = 0; i_sp < nb_ns; i_sp++){ 
//         Mw[i_sp] = mix().speciesMw(i_sp);
//     }
// 
//     int nb_gsi_reactions = mix().nCatalyticReactions();
//     BOOST_CHECK_EQUAL(nb_gsi_reactions, 3); 
// 
//     double gamma_coefficient_NN = 0.01;
//     double gamma_coefficient_OO = 0.02;
//     double gamma_coefficient_NO = 0.03;
//     double gamma_coefficient_ON = 0.04;
//     BOOST_CHECK_CLOSE(gamma_coefficient_NN, 0.01, tol);
//     BOOST_CHECK_CLOSE(gamma_coefficient_OO, 0.02, tol);
//     BOOST_CHECK_CLOSE(gamma_coefficient_NO, 0.03, tol);
//     BOOST_CHECK_CLOSE(gamma_coefficient_ON, 0.04, tol);
// 
//     int nb_pressures = 10;
//     int nb_temperatures = 10;
//     std::vector<double> v_pressure( nb_pressures );
//     std::vector<double> v_temperature( nb_temperatures );
//     get_pressure_temperature_vectors( v_pressure, v_temperature );
// 
//     for (size_t i_pressure = 0 ; i_pressure < nb_pressures; ++i_pressure) {
//         for (size_t i_temperature = 0 ; i_temperature < nb_temperatures; ++i_temperature) {
// 
//             equilibrateMixture(v_temperature[i_temperature], v_pressure[i_pressure]);
//             rhoi = getPartialDensities();
//             mix().setWallState(&rhoi[0],&v_temperature[i_temperature]);
//             mix().netGSIProductionRates( &p_production_rates[0] );
//             sum = sum_production_rates( p_production_rates );
// 
// //            BOOST_CHECK_CLOSE(p_production_rates[1], compute_production_rate_for_O_O2_recombination(gamma_coefficient, Mw, rhoi, v_temperature[i_temperature]), tol );
// //            BOOST_CHECK_CLOSE(sum, 0.0, tol);
// 
//         }
//     }
// 
// }


BOOST_AUTO_TEST_CASE(gsi_reactions_air5_gamma){

    

}

BOOST_AUTO_TEST_SUITE_END()
