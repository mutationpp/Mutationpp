#!/usr/bin/env python

import numpy as np
import mutationpp as mpp
import pytest

#Equilibrium Mixture
mixture = mpp.Mixture("air_11")
mixture.equilibrate(300., 1000.)

# MultiTemperature Mixture
myMixtureOptions = mpp.MixtureOptions("air_5")
myMixtureOptions.setStateModel("ChemNonEqTTv")
mixtureNonEq = mpp.Mixture(myMixtureOptions)
rhoi = [1.0]*5
T = [10000., 300.]
mixtureNonEq.setState(rhoi, T, 1)

# def test_averageDiffusionCoeffs():
# Todo: write proper test

def test_averageHeavyCollisionFreq():
    assert mixture.averageHeavyCollisionFreq() == pytest.approx(4.0906953e7)


#
# def test_averageHeavyThermalSpeed():
# Todo: write proper test

# def test_backwardRateCoefficients():
# Todo: write proper test

def test_density():
    assert mixture.density() == pytest.approx(0.0115663, abs=1e-7)


def test_dRhodP():
    assert mixture.dRhodP() == pytest.approx(1.15663e-05, abs=1e-7)


# def test_diffusionMatrix():
# Todo: write proper test

# def test_dXidT():
# Todo: write proper test
#   a = np.array([1.])
#   mixture.dXidT(a)

def test_electronHeavyCollisionFreq():
    assert mixture.electronHeavyCollisionFreq() == pytest.approx(0.)


# def test_electronMeanFreePath():
# Todo: write proper test

def test_electronThermalConductivity():
    assert mixture.electronThermalConductivity() == pytest.approx(0.)


# def test_electronThermalSpeed():
# Todo: write proper test

# def test_elementPotentialsWrapper():
# Todo: write proper test

# def test_equilDiffFluxFacsP():
# Todo: write proper test

def test_equilibriumSoundSpeed():
    assert mixture.equilibriumSoundSpeed() == pytest.approx(347.765, abs=1e-3)


def test_equilibriumThermalConductivity():
    assert mixture.equilibriumThermalConductivity() == pytest.approx(0.0284695, abs=1e-7)


# def test_forwardRateCoefficients():
# Todo: write proper test

def test_frozenSoundSpeed():
    assert mixture.frozenSoundSpeed() == pytest.approx(347.765, abs=1e-3)


def test_getBField():
    mixture.setBField(1.)
    assert mixture.getBField() == 1.


def test_hasElectrons():
    assert mixture.hasElectrons()


def test_heavyThermalConductivity():
    assert mixture.heavyThermalConductivity() == pytest.approx(0.0209103, abs=1e-7)


def test_internalThermalConductivity():
    assert mixture.internalThermalConductivity(mixture.T()) == pytest.approx(0.00755904, abs=1e-7)


def test_mixtureEnergyMass():
    assert mixture.mixtureEnergyMass() == pytest.approx(-84588.1, abs=1e-3)


def test_mixtureEnergyMole():
    assert mixture.mixtureEnergyMole() == pytest.approx(-2440.39, abs=1e-2)


def test_mixtureEquilibriumCvMass():
    assert mixture.mixtureEquilibriumCvMass() == pytest.approx(722.589, abs=1e-3)


def test_mixtureEquilibriumCpMass():
    assert mixture.mixtureEquilibriumCpMass() == pytest.approx(1010.78, abs=1e-2)


def test_mixtureEquilibriumCpMole():
    assert mixture.mixtureEquilibriumCpMole() == pytest.approx(29.1614, abs=1e-4)


def test_mixtureEquilibriumGamma():
    assert mixture.mixtureEquilibriumGamma() == pytest.approx(1.39883, abs=1e-5)


def test_mixtureFrozenCpMass():
    assert mixture.mixtureFrozenCpMass() == pytest.approx(1010.78, abs=1e-2)


def test_mixtureFrozenCpMole():
    assert mixture.mixtureFrozenCpMole() == pytest.approx(29.1614, abs=1e-4)


def test_mixtureFrozenCvMole():
    assert mixture.mixtureFrozenCvMole() == pytest.approx(20.8469, abs=1e-4)


def test_mixtureFrozenCvMass():
    assert mixture.mixtureFrozenCvMass() == pytest.approx(722.589, abs=1e-3)


def test_mixtureFrozenGamma():
    assert mixture.mixtureFrozenGamma() == pytest.approx(1.39883, abs=1e-5)


def test_mixtureHMass():
    assert mixture.mixtureHMass() == pytest.approx(1869.87, abs=1e-2)


def test_mixtureHMole():
    assert mixture.mixtureHMole() == pytest.approx(53.9465, abs=1e-4)


def test_mixtureMw():
    assert mixture.mixtureMw() == pytest.approx(0.0288503, abs=1e-7)

def test_mixtureSMass():
    assert mixture.mixtureSMass() == pytest.approx(8218.83, abs=1e-2)


def test_meanFreePath():
    assert mixture.meanFreePath() == pytest.approx(1.1470317911509538e-05)

def test_nAtoms():
    assert mixture.nAtoms() == 4


def test_nCondensed():
    assert mixture.nCondensed() == 0


def test_nElements():
    assert mixture.nElements() == 3


def test_nEnergyEqns():
    assert mixture.nEnergyEqns() == 1


# def test_nEquilibriumSteps():
# Todo: write proper test

# def test_nEquilibriumNewtons():
# Todo: write proper test

def test_nGas():
    assert mixture.nGas() == 11


def test_nHeavy():
    assert mixture.nHeavy() == 10


def test_nMassEqns():
    assert mixture.nMassEqns() == 11


def test_nMolecules():
    assert mixture.nMolecules() == 6


def test_nPhases():
    assert mixture.nPhases() == 1


def test_nSpecies():
    assert mixture.nSpecies() == 11


def test_numberDensity():
    assert mixture.numberDensity() == pytest.approx(2.414321232055165e+23)


# def test_omega11ii():
# Todo: write proper test

# def test_omega22ii():
# Todo: write proper test

# def test_phaseMoles():
# Todo: write proper test

def test_P():
    assert mixture.P() == 1000.


def test_reactiveThermalConductivity():
    assert mixture.reactiveThermalConductivity() == pytest.approx(-5.74754e-08)


def test_soretThermalConductivity():
    assert mixture.soretThermalConductivity() == pytest.approx(3.01407e-07)


# def test_stefanMaxwell():
# Todo: write proper test

# def test_speciesCpOverR():
# Todo: write proper test

# def test_speciesHOverRT():
# Todo: write proper test

# def test_speciesMw():
# Todo: write proper test

# def test_speciesSOverR():
# Todo: write proper test

# def test_speciesGOverRT():
# Todo: write proper test

def test_T():
    assert mixture.T() == 300.
    assert mixtureNonEq.T() == 10000.

def test_Tr():
    assert mixtureNonEq.Tr() == 10000.

def test_Tv():
    assert mixtureNonEq.Tv() == 300.

def test_Te():
    assert mixtureNonEq.Te() == 300.

# def test_thermalDiffusionRatios():
# Todo: write proper test

def test_viscosity():
    assert mixture.viscosity() == pytest.approx(1.9426e-5, abs=1e-8)


def test_X():
    # Todo: write proper test
    assert sum(mixture.X()) == 1.


def test_Y():
    # Todo: write proper test
    assert sum(mixture.Y()) == 1.
    
def test_convert_xe_to_ye():
    np.testing.assert_array_equal(
        mixture.convert_xe_to_ye([0, 0, 1]), np.array([0, 0, 1])
    )

def test_convert_ye_to_xe():
    np.testing.assert_array_equal(
        mixture.convert_ye_to_xe([0, 0, 1]), np.array([0, 0, 1])
    )   
    
def test_speciesHOverRT():
    expected = [
        1.54166667e-02,
        7.54484667e02,
        6.28975420e02,
        3.97181733e02,
        6.05194961e02,
        4.69649788e02,
        1.89420120e02,
        9.99331729e01,
        3.65398913e01,
        2.15923313e-02,
        2.17600146e-02,
    ]
    assert np.allclose(mixture.speciesHOverRT(), expected)


def test_speciesHOverRT():
    expected = [
        2.36710963e-02,
        7.51986575e02,
        6.26894106e02,
        3.95873824e02,
        6.03195982e02,
        4.68101138e02,
        1.88799123e02,
        9.96094746e01,
        3.64301530e01,
        3.31535710e-02,
        3.34133823e-02,
    ]
    assert np.allclose(mixture.speciesHOverRT(301), expected)    
    
def test_getEnergiesMass():
    
    expected = mixture.getEnthalpiesMass() -  mixture.T()*mixture.RU()/mixture.speciesMw()
    
    assert mixture.getEnergiesMass() == pytest.approx(expected)    
    
    
def test_getEnthalpiesMass():
    
    expected = [
        7.0098139454e+07,
        1.3436541662e+08,
        9.8062005405e+07,
        3.3017452040e+07,
        5.3888226083e+07,
        3.6610342011e+07,
        3.3732317707e+07,
        1.5579800144e+07,
        3.0374812162e+06,
        1.9226029990e+03,
        1.6962169231e+03
    ]
    
    assert mixture.getEnthalpiesMass() == pytest.approx(expected)  
     
def test_getCpsMass():
    
    expected = [
        3.7890886192e+07,
        1.5192147894e+03,
        1.2992294338e+03,
        9.7019017139e+02,
        1.0395135823e+03,
        9.1103684272e+02,
        1.4840168399e+03,
        1.2991848864e+03,
        9.7219901895e+02,
        1.0392575435e+03,
        9.1700386827e+02
    ]
    
    assert mixture.getCpsMass() == pytest.approx(expected)  
    
    
def test_getGibbsMass():
    
    expected = [
        -3.2412103851e+10,
        1.3011972425e+08,
        9.4435909312e+07,
        3.0651845001e+07,
        5.1361201327e+07,
        3.4318319525e+07,
        2.9626009866e+07,
        1.1827394569e+07,
        5.4328552314e+05,
        -2.4612058750e+06,
        -2.2814767433e+06
    ]
    
    assert mixture.getGibbsMass() == pytest.approx(expected)      
        
def test_getCvsMass():
    
    expected = mixture.getCpsMass()- mixture.RU()/mixture.speciesMw()
    assert mixture.getCvsMass() == pytest.approx(expected)              

def test_mixtureSMass():
    assert mixture.mixtureSMass() == pytest.approx(8218.83)  
    
def test_NA():
    assert mixture.NA() == pytest.approx(6.0221415E23)       
    
def test_KB():
    assert mixture.KB() == pytest.approx(1.3806503E-23)  
    
def test_RU():
    assert mixture.RU() == pytest.approx(mixture.KB()* mixture.NA()) 
    
def test_HP():
    assert mixture.HP() == pytest.approx(6.626068E-34)      
    
def test_C0():
    assert mixture.C0() == pytest.approx(299792458.0)     
    
def test_ONEATM():
    assert mixture.ONEATM() == pytest.approx(101325.0)                               
