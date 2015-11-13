#!/usr/bin/env python

import pyMutation as mpp
import unittest
import numpy as np

# Testing

# class TestMixtureOptions(unittest.TestCase):
#   def setUp(self):
#       self.options = mpp.MixtureOptions("air11")

#   def test_getSpeciesDescriptor(self):
#       self.assertEqual(
#           self.options.getSpeciesDescriptor(),
#           "e- N N+ O O+ NO N2 N2+ O2 O2+ NO+")


class TestMixture(unittest.TestCase):
    # Todo: change the test suite to avoid some errors (see with JB)

    def setUp(self):
        self.mixture = mpp.Mixture("air11")
        self.mixture.equilibrate(300., 1000.)

    # def test_averageDiffusionCoeffs(self):
    # Todo: write proper test

    def test_averageHeavyCollisionFreq(self):
        self.assertAlmostEqual(self.mixture.averageHeavyCollisionFreq(), 0.)

    def test_averageHeavyThermalSpeed(self):
        self.assertAlmostEqual(self.mixture.averageHeavyThermalSpeed(), 0.)

    # def test_backwardRateCoefficients(self):
    # Todo: write proper test

    def test_density(self):
        self.assertAlmostEqual(self.mixture.density(), 0.)

    def test_dRhodP(self):
        self.assertAlmostEqual(self.mixture.dRhodP(), 0.)

    # def test_diffusionMatrix(self):
    # Todo: write proper test

    # def test_dXidT(self):
    # Todo: write proper test
    #   a = np.array([1.])
    #   mixture.dXidT(a)

    def test_electronHeavyCollisionFreq(self):
        self.assertAlmostEqual(self.mixture.electronHeavyCollisionFreq(), 0.)

    def test_electronMeanFreePath(self):
        self.assertAlmostEqual(self.mixture.electronMeanFreePath(), 0.)

    def test_electronThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.electronThermalConductivity(), 0.)

    def test_electronThermalSpeed(self):
        self.assertAlmostEqual(self.mixture.electronThermalSpeed(), 0.)

    # def test_elementPotentialsWrapper(self):
    # Todo: write proper test

    # def test_equilDiffFluxFacsP(self):
    # Todo: write proper test

    def test_equilibriumSoundSpeed(self):
        self.assertAlmostEqual(self.mixture.equilibriumSoundSpeed(), 0.)

    def test_equilibriumThermalConductivity(self):
        self.assertAlmostEqual(
            self.mixture.equilibriumThermalConductivity(), 0.)

    # def test_forwardRateCoefficients(self):
    # Todo: write proper test

    def test_frozenSoundSpeed(self):
        self.assertAlmostEqual(self.mixture.frozenSoundSpeed(), 0.)

    def test_getBField(self):
        self.mixture.setBField(1.)
        self.assertAlmostEqual(self.mixture.getBField(), 1.)

    def test_hasElectrons(self):
        self.assertTrue(self.mixture.hasElectrons())

    def test_heavyThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.heavyThermalConductivity(), 0.)

    def test_internalThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.internalThermalConductivity(), 0.)

    def test_mixtureEnergyMass(self):
        self.assertAlmostEqual(self.mixture.mixtureEnergyMass(), 0.)

    def test_mixtureEnergyMole(self):
        self.assertAlmostEqual(self.mixture.mixtureEnergyMole(), 0.)

    def test_mixtureEquilibriumCvMass(self):
        self.assertAlmostEqual(self.mixture.mixtureEquilibriumCvMass(), 0.)

    def test_mixtureEquilibriumCpMass(self):
        self.assertAlmostEqual(
            self.mixture.mixtureEquilibriumCpMass(), 1131.8, 1)

    def test_mixtureEquilibriumCpMole(self):
        self.assertAlmostEqual(
            self.mixture.mixtureEquilibriumCpMole(), 18.185, 4)

    def test_mixtureEquilibriumGamma(self):
        self.assertAlmostEqual(self.mixture.mixtureEquilibriumGamma(), 0.)

    def test_mixtureFrozenCpMass(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCpMass(), 0.)

    def test_mixtureFrozenCpMole(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCpMole(), 0.)

    def test_mixtureFrozenCvMole(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCvMole(), 0.)

    def test_mixtureFrozenCvMass(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenCvMass(), 0.)

    def test_mixtureFrozenGamma(self):
        self.assertAlmostEqual(self.mixture.mixtureFrozenGamma(), 0.)

    def test_mixtureHMass(self):
        self.assertAlmostEqual(self.mixture.mixtureHMass(), 0.)

    def test_mixtureHMole(self):
        self.assertAlmostEqual(self.mixture.mixtureHMole(), 0.)

    def test_mixtureMw(self):
        self.assertAlmostEqual(self.mixture.mixtureMw(), 0.)

    def test_mixtureSMass(self):
        self.assertAlmostEqual(self.mixture.mixtureSMass(), 0.)

    def test_mixtureSMole(self):
        self.assertAlmostEqual(self.mixture.mixtureSMole(), 0.)

    def test_meanFreePath(self):
        self.assertAlmostEqual(self.mixture.meanFreePath(), 0.)

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
        self.assertAlmostEqual(self.mixture.numberDensity(), 0.)

    # def test_omega11ii(self):
    # Todo: write proper test

    # def test_omega22ii(self):
    # Todo: write proper test

    # def test_phaseMoles(self):
    # Todo: write proper test

    def test_P(self):
        self.assertAlmostEqual(self.mixture.P(), 1000.)

    def test_reactiveThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.reactiveThermalConductivity(), 0.)

    def test_soretThermalConductivity(self):
        self.assertAlmostEqual(self.mixture.soretThermalConductivity(), 0.)

    def test_sigma(self):
        self.assertAlmostEqual(self.mixture.sigma(), 0.)

    def test_sigmaParallel(self):
        self.assertAlmostEqual(self.mixture.sigmaParallel(), 0.)

    def test_sigmaPerpendicular(self):
        self.assertAlmostEqual(self.mixture.sigmaPerpendicular(), 0.)

    def test_sigmaTransverse(self):
        self.assertAlmostEqual(self.mixture.sigmaTransverse(), 0.)

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
        self.assertAlmostEqual(self.mixture.X(), 0.)

    def test_Y(self):
        # Todo: write proper test
        self.assertAlmostEqual(self.mixture.Y(), 0.)

if __name__ == "__main__":

    myMixtureOptions = mpp.MixtureOptions("air11")
    print(myMixtureOptions.getSource())
    # print(myMixtureOptions.getLoadTransport())

    mixture1 = mpp.Mixture("air5")
    mixture1.equilibrate(300., 1000.)
    # D = np.copy(mixture1.diffusionMatrix())
    # x = PrettyTable(D.dtype.names)
    # for row in D:
    #     x.add_row(D)
    # Change some column alignments; default was 'c'
    # x.align['column_one'] = 'r'
    # x.align['col_two'] = 'r'
    # x.align['column_3'] = 'l'
    # print(x)
    a = np.zeros(5)
    mixture1.averageDiffusionCoeffs(a)
    print(a)
    print(mixture1.diffusionMatrix())

    # mixture2 = mpp.Mixture("air11")
    # mixture2.equilibrate(300.,1000.)
    # a = np.copy(mixture1.X())
    # b = np.copy(mixture2.X())
    # print a - b
    # print mixture1.nSpecies()

    # unittest.main()
