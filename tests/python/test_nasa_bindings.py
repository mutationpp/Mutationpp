#!/usr/bin/env python

"""Tests for the NASA polynomial accessors, which expose the fitted
coefficients a symbolic code generator needs. Mutation++ itself only
evaluates the polynomials.
"""

import mutationpp as mpp
import pytest


def make_mixture(database):
    options = mpp.MixtureOptions("air_5")
    options.setStateModel("ChemNonEqTTv")
    options.setThermodynamicDatabase(database)
    return mpp.Mixture(options)


nasa9 = make_mixture("NASA-9")
nasa7 = make_mixture("NASA-7")


@pytest.mark.parametrize("mixture,num_coefficients", [(nasa9, 9), (nasa7, 7)])
def test_shape_of_every_species_polynomial(mixture, num_coefficients):
    for species_index in range(mixture.num_species):
        polynomial = mixture.species_nasa_polynomial(species_index)
        coefficients = polynomial["coefficients"]
        bounds = polynomial["t_bounds"]
        assert len(bounds) == len(coefficients) + 1
        for interval in coefficients:
            assert len(interval) == num_coefficients


@pytest.mark.parametrize("mixture", [nasa9, nasa7])
def test_temperature_bounds_ascend(mixture):
    for species_index in range(mixture.num_species):
        bounds = mixture.species_nasa_polynomial(species_index)["t_bounds"]
        assert bounds == sorted(bounds)
        assert bounds[0] > 0.0


@pytest.mark.parametrize("mixture", [nasa9, nasa7])
def test_coefficients_reproduce_the_evaluated_specific_heat(mixture):
    # cp/R for NASA-7 is a polynomial in T; for NASA-9 it carries two
    # inverse-power terms as well. Rebuilding it from the coefficients
    # must reproduce what Mutation++ evaluates.
    temperature = 3000.0
    evaluated = mixture.speciesCpOverR(temperature)
    for species_index in range(mixture.num_species):
        polynomial = mixture.species_nasa_polynomial(species_index)
        bounds = polynomial["t_bounds"]
        interval = min(
            sum(temperature > bound for bound in bounds[1:-1]),
            len(polynomial["coefficients"]) - 1,
        )
        c = polynomial["coefficients"][interval]
        if len(c) == 7:
            rebuilt = (c[0] + c[1] * temperature + c[2] * temperature**2
                       + c[3] * temperature**3 + c[4] * temperature**4)
        else:
            rebuilt = (c[0] / temperature**2 + c[1] / temperature + c[2]
                       + c[3] * temperature + c[4] * temperature**2
                       + c[5] * temperature**3 + c[6] * temperature**4)
        assert rebuilt == pytest.approx(evaluated[species_index], rel=1e-10)


def test_rrho_database_has_no_polynomials():
    mixture = make_mixture("RRHO")
    with pytest.raises(RuntimeError, match="NASA"):
        mixture.species_nasa_polynomial(0)
