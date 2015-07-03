
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE catalysis tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "MixtureFixtures.h"

BOOST_FIXTURE_TEST_SUITE(catalysis_wall_solver, MixtureFromCommandLine)

BOOST_AUTO_TEST_CASE(O2_rhois_at_wall)
{
    const double tol = 1.0e3*std::numeric_limits<double>::epsilon();

    const int nb_ns = mix().nSpecies();
    const int nt = mix().nEnergyEqns();
    std::vector<double> rhoi(nb_ns);
    std::vector<double> Mw(nb_ns);
    std::vector<double> p_production_rates(nb_ns);
    double sum = 0.0;
       
    for (int i_sp = 0; i_sp < nb_ns; i_sp++){ 
        Mw[i_sp] = mix().speciesMw(i_sp);
    }

//    const size_t iO   = mix().speciesIndex("O");
//    const size_t iO2  = mix().speciesIndex("O2");

    const double T_free_stream = 1000.E0;
    const double P_free_stream = 100.E0;
    const double T_wall        = 300.E0;

    equilibrateMixture(T_free_stream, P_free_stream);
    // Get the species densities
    double rho = mix().density();
    for (int i = 0; i < nb_ns; ++i)
        rhoi[i] = mix().Y()[i] * rho;

    double discreteLengthforDerivative = 1.E-2;

    std::cout << "First rhois are = " << rhoi[0] << "      " << rhoi[1] << std::endl;
    //mix().solveWallMassBalance(&rhoi[0], &rhoi[0], &T_wall, discreteLengthforDerivative);
    mix().solveSurfaceBalance();
    std::cout << "Then rhoi becomes = " << rhoi[0] <<"      " << rhoi[1] << std::endl;

//       BOOST_CHECK_CLOSE(p_production_rates[0], compute_production_rate_for_O_O2_recombination(gamma_coefficient, Mw, rhoi, v_temperature[i_temperature]), tol );
//       BOOST_CHECK_CLOSE(sum, 0.0, tol);

}

BOOST_AUTO_TEST_SUITE_END()


