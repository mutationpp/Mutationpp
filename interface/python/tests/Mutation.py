#!/usr/bin/env python

import mutationpp_python as mpp
import unittest


class TestMixture(unittest.TestCase):
    # Todo: change the test suite to avoid some errors (see with JB)

    def setUp(self):
        self.mixture = mpp.Mixture("air_11")
        self.mixture.equilibrate(300., 1000.)

    # def test_averageDiffusionCoeffs(self):
    # Todo: write proper test

    def test_averageHeavyCollisionFreq(self):
        self.assertAlmostEqual(self.mixture.averageHeavyCollisionFreq(), 0.)

    # def test_averageHeavyThermalSpeed(self):
    # Todo: write proper test

    # def test_backwardRateCoefficients(self):
    # Todo: write proper test

    def test_density(self):
        self.assertAlmostEqual(self.mixture.density(), 0.0115663)

    def test_dRhodP(self):
        self.assertAlmostEqual(self.mixture.dRhodP(), 1.15663e-05)

    # def test_diffusionMatrix(self):
    # Todo: write proper test

    # def test_dXidT(self):
    # Todo: write proper test
    #   a = np.array([1.])
    #   mixture.dXidT(a)

    def test_electronHeavyCollisionFreq(self):
        self.assertAlmostEqual(self.mixture.electronHeavyCollisionFreq(), 0.)

    # def test_electronMeanFreePath(self):
    # Todo: write proper test

    def test_electronThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.electronThermalConductivity(), 0.)

    # def test_electronThermalSpeed(self):
    # Todo: write proper test

    # def test_elementPotentialsWrapper(self):
    # Todo: write proper test

    # def test_equilDiffFluxFacsP(self):
    # Todo: write proper test

    def test_equilibriumSoundSpeed(self):
        self.assertAlmostEqual(self.mixture.equilibriumSoundSpeed(), 347.765, 3)

    def test_equilibriumThermalConductivity(self):
        self.assertAlmostEqual(
            self.mixture.equilibriumThermalConductivity(), 0.0284695, 7)

    # def test_forwardRateCoefficients(self):
    # Todo: write proper test

    def test_frozenSoundSpeed(self):
        self.assertAlmostEqual(self.mixture.frozenSoundSpeed(), 347.765, 3)

    def test_getBField(self):
        self.mixture.setBField(1.)
        self.assertAlmostEqual(self.mixture.getBField(), 1.)

    def test_hasElectrons(self):
        self.assertTrue(self.mixture.hasElectrons())

    def test_heavyThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.heavyThermalConductivity(), 0.0209103, 7)

    def test_internalThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.internalThermalConductivity(self.mixture.T()), 0.00755904, 8)

    def test_mixtureEnergyMass(self):
        self.assertAlmostEqual(self.mixture.mixtureEnergyMass(), -84588.1, 2)

    def test_mixtureEnergyMole(self):
        self.assertAlmostEqual(self.mixture.mixtureEnergyMole(), -2440.39, 2)

    def test_mixtureEquilibriumCvMass(self):
        self.assertAlmostEqual(self.mixture.mixtureEquilibriumCvMass(), 722.589, 3)

    def test_mixtureEquilibriumCpMass(self):
        self.assertAlmostEqual(
            self.mixture.mixtureEquilibriumCpMass(), 1010.78, 2)

    def test_mixtureEquilibriumCpMole(self):
        self.assertAlmostEqual(
            self.mixture.mixtureEquilibriumCpMole(), 29.1614, 4)

    def test_mixtureEquilibriumGamma(self):
        self.assertAlmostEqual(self.mixture.mixtureEquilibriumGamma(), 1.39883, 4)

    def test_mixtureFrozenCpMass(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCpMass(), 1010.78, 2)

    def test_mixtureFrozenCpMole(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCpMole(), 29.1614, 4)

    def test_mixtureFrozenCvMole(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCvMole(), 20.8469, 4)

    def test_mixtureFrozenCvMass(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCvMass(), 722.589, 3)

    def test_mixtureFrozenGamma(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenGamma(), 1.39883, 4)

    def test_mixtureHMass(self):
        self.assertAlmostEqual(self.mixture.mixtureHMass(), 1869.87, 2)

    def test_mixtureHMole(self):
        self.assertAlmostEqual(self.mixture.mixtureHMole(), 53.9465, 4)

    def test_mixtureMw(self):
        self.assertAlmostEqual(self.mixture.mixtureMw(), 0.0288503)

    # def test_meanFreePath(self):
    #     Todo: write proper test

    def test_nAtoms(self):
        self.assertEqual(self.mixture.nAtoms(), 4)

    def test_nCondensed(self):
        self.assertEqual(self.mixture.nCondensed(), 0)

    def test_nElements(self):
        self.assertEqual(self.mixture.nElements(), 3)

    def test_nEnergyEqns(self):
        self.assertEqual(self.mixture.nEnergyEqns(), 1)

    # def test_nEquilibriumSteps(self):
    # Todo: write proper test

    # def test_nEquilibriumNewtons(self):
    # Todo: write proper test

    def test_nGas(self):
        self.assertEqual(self.mixture.nGas(), 11)

    def test_nHeavy(self):
        self.assertEqual(self.mixture.nHeavy(), 10)

    def test_nMassEqns(self):
        self.assertEqual(self.mixture.nMassEqns(), 11)

    def test_nMolecules(self):
        self.assertEqual(self.mixture.nMolecules(), 6)

    def test_nPhases(self):
        self.assertEqual(self.mixture.nPhases(), 1)

    def test_nSpecies(self):
        self.assertEqual(self.mixture.nSpecies(), 11)

    def test_numberDensity(self):
        self.assertAlmostEqual(self.mixture.numberDensity(), 2.414321232055165e+23)

    # def test_omega11ii(self):
    # Todo: write proper test

    # def test_omega22ii(self):
    # Todo: write proper test

    # def test_phaseMoles(self):
    # Todo: write proper test

    def test_P(self):
        self.assertAlmostEqual(self.mixture.P(), 1000.)

    def test_reactiveThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.reactiveThermalConductivity(), -5.74754e-08)

    def test_soretThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.soretThermalConductivity(), 3.01407e-07)

    # def test_stefanMaxwell(self):
    # Todo: write proper test

    # def test_speciesCpOverR(self):
    # Todo: write proper test

    # def test_speciesHOverRT(self):
    # Todo: write proper test

    # def test_speciesMw(self):
    # Todo: write proper test

    # def test_speciesSOverR(self):
    # Todo: write proper test

    # def test_speciesGOverRT(self):
    # Todo: write proper test

    def test_T(self):
        self.assertAlmostEqual(self.mixture.T(), 300.)

    # def test_thermalDiffusionRatios(self):
    # Todo: write proper test

    def test_viscosity(self):
        self.assertAlmostEqual(self.mixture.viscosity(), 1.9426e-5)

    def test_X(self):
        # Todo: write proper test
        self.assertAlmostEqual(sum(self.mixture.X()), 1.)

    def test_Y(self):
        # Todo: write proper test
        self.assertAlmostEqual(sum(self.mixture.Y()), 1.)


if __name__ == "__main__":
    unittest.main()
