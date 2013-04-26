/**
 * @file ThermoTests.h
 *
 * Implements common tests for thermodynamic data which can be applied to any
 * mixture.  Each test can be turned on by definining the appropriate 
 * preprocessor variable before including this file.
 */
#include <stdio.h>
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
#define TESTE EXPAND_AND_QUOTE(TEST_DATA_DIRECTORY)


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
        TESTE "/" TEST_FIXTURE_DATA_DIRECTORY "/" #__FUNC__ ".dat", __FUNC__, __NVALS__ );\
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


//  14: Cp        [J/kg-K]    frozen specific heat at constant pressure
// double Mutation::Thermodynamics::Thermodynamics::mixtureFrozenCpMass() const 	
#ifdef TEST_MIXTURE_CP
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_CP, 1, 
    values[0] = mix->mixtureFrozenCpMass();
)
#endif



//  8 : Cp_eq     [J/kg-K]    equilibrium specific heat at constant pressure
// double Mutation::Thermodynamics::Thermodynamics::mixtureEquilibriumCpMass()
#ifdef TEST_MIXTURE_CP_EQ
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_CP_EQ, 1, 
    values[1] = mix->mixtureEquilibriumCpMass();
)
#endif



//  8 : Cp_eq     [J/kg-K]    equilibrium specific heat at constant pressure
// double Mutation::Thermodynamics::Thermodynamics::mixtureEquilibriumCpMass()
#ifdef TEST_MIXTURE_CV
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_CV, 1, 
    values[1] = mix->mixtureFrozenCvMass();
)
#endif


//  8 : Cp_eq     [J/kg-K]    equilibrium specific heat at constant pressure
// double Mutation::Thermodynamics::Thermodynamics::mixtureEquilibriumCpMass()
#ifdef TEST_MIXTURE_CV_EQ
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_CV_EQ, 1, 
    values[1] = mix->mixtureEquilibriumCvMass();
)
#endif




// "Mixture Enthalpy"  9 : H         [J/kg]      mixture enthalpy
// double Mutation::Thermodynamics::Thermodynamics::mixtureHMass() const
#ifdef TEST_MIXTURE_H
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_H, 1, 
    values[0] = mix->mixtureHMass();
)
#endif



// "Mixture Entropy"  10: S         [J/kg-K]    entropy
// double Mutation::Thermodynamics::Thermodynamics::mixtureSMass() const
#ifdef TEST_MIXTURE_S
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_S, 1, 
    values[0] = mix->mixtureSMass();
)
#endif


// "Mixture Energy" 29: e         [J/kg]      mixture energy
//double Mutation::Thermodynamics::Thermodynamics::mixtureEnergyMass() const
#ifdef TEST_MIXTURE_E
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_E, 1, 
 values[0] = mix->mixtureEnergyMass();
)
#endif


// "Rotational Enthalpy" 24: Hr        [J/kg]      rotational enthalpy
// double Mutation::Thermodynamics::ParticleRRHO::rotationalTemperature() const
#ifdef TEST_MIXTURE_HR
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_HR, 1, 
    values[0] = mix->rotationalTemperature();
)
#endif



// "Formation Enthalpy " 27: Hf        [J/kg]      formation enthalpy
// double Mutation::Thermodynamics::ParticleRRHO::formationEnthalpy() 	const
#ifdef TEST_MIXTURE_HF
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_HF, 1, 
    values[0] = mix->formationEnthalpy();
)
#endif


// "Electronic Enthalpy" 26: Hel       [J/kg]      electronic enthalpy
// const std::pair<int, double>& Mutation::Thermodynamics::ParticleRRHO::electronicEnergy(const int i)	const 
#ifdef TEST_MIXTURE_HEL
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_HEL, 1, 
    values[0] = mix->electronicEnergy();
)
#endif


// "Vibration Enthalpy" 25: Hv        [J/kg]      vibrational enthalpy
// double Mutation::Thermodynamics::ParticleRRHO::vibrationalEnergy(const int i) const 
#ifdef TEST_MIXTURE_HV
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_HV, 1, 
    values[0] = mix->vibrationalEnergy();
)
#endif


// "Mixture Frozen Sound Speed" 37: a_f       [m/s]       frozen speed of sound
// double Mutation::Thermodynamics::Thermodynamics::frozenSoundSpeed() const
#ifdef TEST_MIXTURE_A_F
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_A_F, 1, 
    values[0] = mix->frozenSoundSpeed();
)
#endif


// "Mixture Equilibrium Sound Speed" 38: a_eq      [m/s]       equilibrium speed of sound
// double Mutation::Thermodynamics::Thermodynamics::equilibriumSoundSpeed() 
#ifdef TEST_MIXTURE_A_EQ
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_A_EQ, 1, 
    values[0] = mix->equilibriumSoundSpeed();
)
#endif

// "Mixture Frozen Specific Heat Ratio" 17: gamma     [-]         frozen ratio of specific heat
// double Mutation::Thermodynamics::Thermodynamics::mixtureFrozenGamma() const
#ifdef TEST_MIXTURE_GAMMA
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_GAMMA, 1, 
    values[0] = mix->mixtureFrozenGamma();
)
#endif


//  "Mixture Equilibrium Specific Heat Ratio" 16: gam_eq    [-]         equilibrium ratio of specific heats
// double Mutation::Thermodynamics::Thermodynamics::mixtureEquilibriumGamma() 
#ifdef TEST_MIXTURE_GAM_EQ
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_GAM_EQ, 1, 
    values[0] = mix->mixtureEquilibriumGamma();
)
#endif


#ifdef TEST_MIXTURE_RHO
TEST_COMPARE_EQUILIBRIUM_VALUES(MIX_RHO, 1, 
    values[0] = mix->density();
)
#endif




#ifdef TEST_SPECIES_X
TEST_COMPARE_EQUILIBRIUM_VALUES(SPE_X,  mix()->nSpecies(),
    for (int i = 0; i < ns; ++i)
        values[i] = mix->X()[i];
)
#endif


#ifdef TEST_SPECIES_Y
TEST_COMPARE_EQUILIBRIUM_VALUES(SPE_Y, mix()->nSpecies(),
    for (int i = 0; i < ns; ++i)
        values[i] = mix->Y()[i];
)
#endif






BOOST_AUTO_TEST_SUITE_END()
/**



// TEST TEMPLATE
//#ifdef TEST_MIXTURE_(NAME)
//TEST_COMPARE_EQUILIBRIUM_VALUES((NAME), (# of VALUES), 
//    values[0] = mix->(function call);
//)
//#endif

// End the test suite
BOOST_AUTO_TEST_SUITE_END()


*/

