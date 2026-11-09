#!/usr/bin/env python

"""Tests for the mode-resolved enthalpy binding, which splits each
species' enthalpy into its translational, rotational, vibrational,
electronic and formation parts.
"""

import mutationpp as mpp
import pytest


options = mpp.MixtureOptions("air_5")
options.setStateModel("ChemNonEqTTv")
mixture = mpp.Mixture(options)
mixture.setState([1e-3] * 5, [9000.0, 4000.0], 1)

species_names = [
    mixture.speciesName(index) for index in range(mixture.nSpecies())
]
modes = [
    "translational", "rotational", "vibrational", "electronic", "formation"
]


def test_every_mode_is_returned():
    enthalpies = mixture.species_enthalpies_over_rt()
    assert set(enthalpies) == set(modes) | {"total"}
    for values in enthalpies.values():
        assert len(values) == mixture.nSpecies()


def test_modes_sum_to_the_total():
    enthalpies = mixture.species_enthalpies_over_rt()
    for index in range(mixture.nSpecies()):
        assert sum(enthalpies[mode][index] for mode in modes) == (
            pytest.approx(enthalpies["total"][index])
        )


def test_total_agrees_with_the_aggregate_accessor():
    enthalpies = mixture.species_enthalpies_over_rt()
    for index, value in enumerate(mixture.speciesHOverRT()):
        assert enthalpies["total"][index] == pytest.approx(value)


def test_atoms_carry_no_rotational_or_vibrational_enthalpy():
    enthalpies = mixture.species_enthalpies_over_rt()
    for species_name in ["N", "O"]:
        index = species_names.index(species_name)
        assert enthalpies["rotational"][index] == pytest.approx(0.0)
        assert enthalpies["vibrational"][index] == pytest.approx(0.0)
        assert enthalpies["electronic"][index] != pytest.approx(0.0)


def test_translational_enthalpy_is_five_halves():
    enthalpies = mixture.species_enthalpies_over_rt()
    for index in range(mixture.nSpecies()):
        assert enthalpies["translational"][index] == pytest.approx(2.5)
