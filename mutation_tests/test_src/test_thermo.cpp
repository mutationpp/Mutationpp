
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Thermo
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "MppTestFixture.h"
using namespace std;
void cp(Mutation::Mixture* const p_mix, double* const values) {
    values[0] = p_mix->mixtureFrozenCpMass();
    values[1] = p_mix->mixtureEquilibriumCpMass();

}

// Air11 RRHO
MPP_TEST_FIXTURE( Air11RRHO, "Air11-RRHO", "air11",
    options.setThermodynamicDatabase("RRHO");
)

BOOST_FIXTURE_TEST_SUITE(Thermodynamics, Air11RRHO)


BOOST_AUTO_TEST_CASE(MixtureCp)
{
    string path = boost::filesystem::current_path().generic_string();

cout << endl << endl<< endl  "Path ==> " << path << endl << endl << endl; 
    compareEquilibriumValues("/home/didi/mutation++/branches/dinesh/dinesh/mutation_tests/test_src/air11/RRHO/cp.dat", cp, 2);
}



BOOST_AUTO_TEST_CASE(SetState)
{
    double T = 1000.0;
    double P = 101325.0;
    
    std::fill(sp1(), sp1()+mix()->nSpecies(), 1.0 / mix()->nSpecies());
    mix()->setStateTPX(&T, &P, sp1());
    
    BOOST_CHECK(mix()->T() == 1000.0);
    BOOST_CHECK(mix()->P() == 101325.0);
    cout << "mix()->nSpecies() : " << mix()->nSpecies() << endl;
    for (int i = 0; i < mix()->nSpecies(); ++i)
        BOOST_CHECK_EQUAL(mix()->X()[i] , 1.0 / mix()->nSpecies());
}

BOOST_AUTO_TEST_SUITE_END()




