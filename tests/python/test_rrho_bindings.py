#!/usr/bin/env python

"""Tests for the ParticleRRHO bindings, which expose the rigid-rotor
harmonic-oscillator parameters a symbolic code generator needs: the
NASA databases cannot be used for two-temperature work because they
zero out the vibrational and electronic specific heats.
"""

import mutationpp as mpp
import pytest


mixture = mpp.Mixture("air_5")
species_names = [
    mixture.speciesName(index) for index in range(mixture.nSpecies())
]


def rrho(species_name):
    return mixture.species_rrho(species_names.index(species_name))


@pytest.mark.parametrize("species_name", ["N", "O", "NO", "N2", "O2"])
def test_every_species_has_electronic_levels(species_name):
    levels = rrho(species_name).electronic_levels
    assert levels
    for degeneracy, characteristic_temperature in levels:
        assert degeneracy >= 1
        assert characteristic_temperature >= 0.0


def test_atoms_have_no_vibrational_mode():
    for species_name in ["N", "O"]:
        assert rrho(species_name).vibrational_temperatures == []
        assert rrho(species_name).linearity == 0


def test_molecules_are_linear_with_one_vibrational_temperature():
    # air5's diatomics, against data/thermo/species.xml.
    expected = {"NO": 2759.293, "N2": 3408.464, "O2": 2276.979}
    for species_name, characteristic_temperature in expected.items():
        particle = rrho(species_name)
        assert particle.linearity == 2
        assert particle.vibrational_temperatures == pytest.approx(
            [characteristic_temperature]
        )


def test_rotational_temperatures_match_the_database():
    expected = {"NO": 2.464, "N2": 2.886, "O2": 2.086}
    for species_name, rotational_temperature in expected.items():
        assert rrho(species_name).rotational_temperature == pytest.approx(
            rotational_temperature
        )


def test_electronic_level_counts_match_the_database():
    assert len(rrho("N").electronic_levels) == 46
    assert len(rrho("O").electronic_levels) == 40


def test_ground_electronic_level_is_at_zero_energy():
    for species_name in species_names:
        degeneracy, energy = rrho(species_name).electronic_levels[0]
        assert energy == pytest.approx(0.0)
        assert degeneracy >= 1


def test_formation_enthalpy_ordering_is_physical():
    # N and O are formed by breaking bonds, so they sit well above the
    # diatomics; N2 and O2 are the reference species at zero.
    assert rrho("N2").formation_enthalpy == pytest.approx(0.0, abs=1.0)
    assert rrho("O2").formation_enthalpy == pytest.approx(0.0, abs=1.0)
    assert rrho("N").formation_enthalpy > rrho("NO").formation_enthalpy > 0.0
