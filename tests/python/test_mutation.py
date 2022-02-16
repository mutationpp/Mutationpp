#!/usr/bin/env python

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
    assert mixture.averageHeavyCollisionFreq() == pytest.approx(0.0)


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


# def test_meanFreePath():
#     Todo: write proper test

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
