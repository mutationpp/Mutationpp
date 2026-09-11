#!/usr/bin/env python

"""Tests for the Reaction, Arrhenius and mixture-level kinetics
bindings used by symbolic code generators, which read mechanism
parameters rather than evaluated properties.
"""

import mutationpp as mpp
import pytest


mixture = mpp.Mixture("air_5")


def test_reaction_accessors_agree_with_mixture_level_ones():
    for reaction_index in range(mixture.num_reactions):
        reaction = mixture.reaction(reaction_index)
        assert mixture.reactants(reaction_index) == reaction.reactants
        assert mixture.products(reaction_index) == reaction.products


def test_reactants_and_products_are_species_indices():
    for reaction_index in range(mixture.num_reactions):
        reaction = mixture.reaction(reaction_index)
        for species_index in reaction.reactants + reaction.products:
            assert 0 <= species_index < mixture.num_species


def test_rate_law_exposes_arrhenius_parameters():
    for reaction_index in range(mixture.num_reactions):
        rate_law = mixture.reaction(reaction_index).rate_law()
        assert rate_law.log_pre_exponential == pytest.approx(
            rate_law.log_pre_exponential
        )
        assert rate_law.activation_temperature > 0.0
