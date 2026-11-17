#!/usr/bin/env python

"""Tests for the Millikan-White bindings, which expose the
vibration-translation relaxation parameters of each vibrator.
"""

import mutationpp as mpp
import pytest


mixture = mpp.Mixture("air_5")
species_names = [
    mixture.speciesName(index) for index in range(mixture.nSpecies())
]
vibrators = ["NO", "N2", "O2"]


def millikan_white(species_name):
    return mixture.millikan_white_data(species_names.index(species_name))


@pytest.mark.parametrize("species_name", vibrators)
def test_one_coefficient_pair_per_heavy_partner(species_name):
    data = millikan_white(species_name)
    assert len(data.a) == mixture.nHeavy()
    assert len(data.b) == mixture.nHeavy()


@pytest.mark.parametrize("species_name", vibrators)
def test_coefficients_are_positive(species_name):
    data = millikan_white(species_name)
    assert all(coefficient > 0.0 for coefficient in data.a)
    assert all(coefficient > 0.0 for coefficient in data.b)


@pytest.mark.parametrize("species_name", vibrators)
def test_species_index_and_molecular_weight_round_trip(species_name):
    data = millikan_white(species_name)
    species_index = species_names.index(species_name)
    assert data.species_index == species_index
    assert data.molecular_weight == pytest.approx(
        mixture.speciesMw(species_index)
    )


def test_database_overrides_the_analytic_defaults():
    # data/transfer/VT.xml gives N2 explicit Millikan-White constants,
    # so they must not equal the 1.16e-3 sqrt(mu) theta^4/3 default.
    nitrogen = millikan_white("N2")
    assert nitrogen.a[species_names.index("N2")] == pytest.approx(221.0)
    assert nitrogen.b[species_names.index("N2")] == pytest.approx(0.0290)
    assert nitrogen.a[species_names.index("O")] == pytest.approx(72.4)


def test_reference_cross_section_is_set_from_the_database():
    # VT.xml declares omegav=3.0e-21 for N2, against the 1.0e-20 default.
    assert millikan_white("N2").reference_cross_section == pytest.approx(
        3.0e-21
    )
