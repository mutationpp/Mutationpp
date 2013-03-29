/**
 * @file ThermoTests.h
 *
 * Implements common tests for thermodynamic data which can be applied to any
 * mixture.  Each test can be turned on by definining the appropriate 
 * preprocessor variable before including this file.
 */
#include <stdio.h>



#define VALUES_TO_COMPARE(__FUNC__,__VALUES__)\
void __FUNC__ (MppTestFixture* const p_fixture, double* const values) {\
    Mutation::Mixture* mix = p_fixture->mix();\
    double* const sp1 = p_fixture->sp1();\
    const int ns = mix->nSpecies();\
    __VALUES__\
}

#define TEST_COMPARE_EQUILIBRIUM_VALUES(__FUNC__,__NVALS__,__VALUES__)\
VALUES_TO_COMPARE( __FUNC__, __VALUES__ )\
BOOST_AUTO_TEST_CASE( compare_equilibrium_##__FUNC__ )\
{\
    compareEquilibriumValues(\
        TEST_DATA_DIRECTORY "/" TEST_FIXTURE_DATA_DIRECTORY "/" #__FUNC__ ".dat", __FUNC__, __NVALS__ );\
}

// Start the test suite
BOOST_FIXTURE_TEST_SUITE(Thermodynamics, TEST_FIXTURE_NAME)

#ifdef TEST_SET_STATE
BOOST_AUTO_TEST_CASE(SetState)
{



    double T = 1000.0;
    double P = 101325.0;
    
    std::fill(sp1(), sp1()+mix()->nSpecies(), 1.0 / mix()->nSpecies());
    mix()->setStateTPX(&T, &P, sp1());
    
    BOOST_CHECK_EQUAL(mix()->T(), 1000.0);
    BOOST_CHECK_EQUAL(mix()->P(), 101325.0);
    
    for (int i = 0; i < mix()->nSpecies(); ++i)
        BOOST_CHECK_EQUAL(mix()->X()[i], 1.0 / mix()->nSpecies());
}
#endif

#ifdef TEST_MIXTURE_CP
TEST_COMPARE_EQUILIBRIUM_VALUES(cp, 2, 
    values[0] = mix->mixtureFrozenCpMass();
    values[1] = mix->mixtureEquilibriumCpMass();
)
#endif

#ifdef TEST_MIXTURE_H
TEST_COMPARE_EQUILIBRIUM_VALUES(h, 1, 
    values[0] = mix->mixtureHMass();
)
#endif

#ifdef TEST_SPECIES_X
TEST_COMPARE_EQUILIBRIUM_VALUES(X, mix()->nSpecies(),
    for (int i = 0; i < ns; ++i)
        values[i] = mix->X()[i];
)
#endif

// TEST TEMPLATE
//#ifdef TEST_MIXTURE_(NAME)
//TEST_COMPARE_EQUILIBRIUM_VALUES((NAME), (# of VALUES), 
//    values[0] = mix->(function call);
//)
//#endif

// End the test suite
BOOST_AUTO_TEST_SUITE_END()


