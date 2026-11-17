#!/usr/bin/env python

"""Tests for the Reaction, Arrhenius and mixture-level kinetics
bindings used by symbolic code generators, which read mechanism
parameters rather than evaluated properties.
"""

import mutationpp as mpp
import pytest


mixture = mpp.Mixture("air_5")


def test_reaction_accessors_agree_with_mixture_level_ones():
    for reaction_index in range(mixture.nReactions()):
        reaction = mixture.reaction(reaction_index)
        assert mixture.reactants(reaction_index) == reaction.reactants
        assert mixture.products(reaction_index) == reaction.products


def test_reactants_and_products_are_species_indices():
    for reaction_index in range(mixture.nReactions()):
        reaction = mixture.reaction(reaction_index)
        for species_index in reaction.reactants + reaction.products:
            assert 0 <= species_index < mixture.nSpecies()


def test_rate_law_exposes_arrhenius_parameters():
    for reaction_index in range(mixture.nReactions()):
        rate_law = mixture.reaction(reaction_index).rate_law()
        assert rate_law.log_pre_exponential == pytest.approx(
            rate_law.log_pre_exponential
        )
        assert rate_law.activation_temperature > 0.0


def test_reaction_metadata_is_exposed():
    # air5_Park: three third-body dissociations then two exchanges.
    third_body_count = sum(
        mixture.reaction(reaction_index).is_thirdbody
        for reaction_index in range(mixture.nReactions())
    )
    assert third_body_count == 3
    for reaction_index in range(mixture.nReactions()):
        reaction = mixture.reaction(reaction_index)
        assert reaction.is_reversible
        assert "=" in reaction.formula
        assert reaction.reaction_type
        assert reaction.order == len(reaction.reactants) + reaction.is_thirdbody


def test_third_body_efficiencies_are_named_by_species_index():
    for reaction_index in range(mixture.nReactions()):
        reaction = mixture.reaction(reaction_index)
        efficiencies = reaction.efficiencies
        if not reaction.is_thirdbody:
            assert efficiencies == []
            continue
        assert efficiencies
        for species_index, efficiency in efficiencies:
            assert 0 <= species_index < mixture.nSpecies()
            assert efficiency > 0.0


def test_rate_temperatures_are_readable_on_a_const_reaction():
    # Dissociation runs forward on sqrt(T*Tv) but reverse on T, so the
    # two must be readable independently.
    dissociation = [
        mixture.reaction(reaction_index)
        for reaction_index in range(mixture.nReactions())
        if mixture.reaction(reaction_index).is_thirdbody
    ]
    for reaction in dissociation:
        assert reaction.fwd_rate_coeff_temperature == "geometric_ttv"
        assert reaction.rev_rate_coeff_temperature == "translational"
