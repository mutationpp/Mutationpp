/**
 * @file RateManager.cpp
 *
 * @brief Implementation of RateManager class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <typeinfo>

#include "RateManager.h"
#include "Reaction.h"
#include "StateModel.h"

namespace Mutation {
    namespace Kinetics {

//==============================================================================
    
// Simple macro to create a temperature selector type
#define TEMPERATURE_SELECTOR(__NAME__,__T__)\
class __NAME__\
{\
public:\
    inline double getT(const Thermodynamics::StateModel* const state) const {\
        return ( __T__ );\
    }\
};

/// Temperature selector which returns the current translational temperature
TEMPERATURE_SELECTOR(TSelector, state->T())

/// Temperature selector which returns the current electron temperature
//TEMPERATURE_SELECTOR(TeSelector, std::min(state->Te(), 10000.0))
TEMPERATURE_SELECTOR(TeSelector, state->Te())

/// Temperature selector which returns the current value of sqrt(T*Tv)
TEMPERATURE_SELECTOR(ParkSelector, std::sqrt(state->T()*state->Tv()))

#undef TEMPERATURE_SELECTOR

/// Arrhenius group evaluated at T
typedef RateLawGroup1T<Arrhenius, TSelector> ArrheniusT;

/// Arrhenius group evaluated at Te
typedef RateLawGroup1T<Arrhenius, TeSelector> ArrheniusTe;

/// Arrhenius group evaluated at sqrt(T*Tv)
typedef RateLawGroup1T<Arrhenius, ParkSelector> ArrheniusPark; 

//==============================================================================

/**
 * Used to define which rate law groups the forward and reverse rate laws are
 * evaluated in for a given ReactionType value.  The default is an Arrhenius
 * rate law with Tf = Tb = T.
 */
template <int Type>
struct RateSelector {
    typedef ArrheniusT ForwardGroup;
    typedef ArrheniusT ReverseGroup;
};

#define SELECT_RATE_LAWS(__TYPE__,__FORWARD__,__REVERSE__)\
template <> struct RateSelector<__TYPE__> {\
    typedef __FORWARD__ ForwardGroup;\
    typedef __REVERSE__ ReverseGroup;\
};

// Default rate law groups for non (kf(T), kb(T)) reaction types
SELECT_RATE_LAWS(ASSOCIATIVE_IONIZATION,     ArrheniusT,    ArrheniusTe)
SELECT_RATE_LAWS(DISSOCIATIVE_RECOMBINATION, ArrheniusTe,   ArrheniusT)
SELECT_RATE_LAWS(ASSOCIATIVE_DETACHMENT,     ArrheniusT,    ArrheniusTe)
SELECT_RATE_LAWS(DISSOCIATIVE_ATTACHMENT,    ArrheniusTe,   ArrheniusT)
SELECT_RATE_LAWS(DISSOCIATION_E,             ArrheniusTe,   ArrheniusTe)
SELECT_RATE_LAWS(RECOMBINATION_E,            ArrheniusTe,   ArrheniusTe)
SELECT_RATE_LAWS(DISSOCIATION_M,             ArrheniusPark, ArrheniusT)
SELECT_RATE_LAWS(RECOMBINATION_M,            ArrheniusT,    ArrheniusPark)
SELECT_RATE_LAWS(IONIZATION_E,               ArrheniusTe,   ArrheniusT)
SELECT_RATE_LAWS(ION_RECOMBINATION_E,        ArrheniusT,    ArrheniusTe)
SELECT_RATE_LAWS(IONIZATION_M,               ArrheniusT,    ArrheniusT)
SELECT_RATE_LAWS(ION_RECOMBINATION_M,        ArrheniusT,    ArrheniusT)
SELECT_RATE_LAWS(ELECTRONIC_ATTACHMENT_M,    ArrheniusTe,   ArrheniusT)
SELECT_RATE_LAWS(ELECTRONIC_DETACHMENT_M,    ArrheniusT,    ArrheniusTe)
SELECT_RATE_LAWS(ELECTRONIC_ATTACHMENT_E,    ArrheniusTe,   ArrheniusTe)
SELECT_RATE_LAWS(ELECTRONIC_DETACHMENT_E,    ArrheniusTe,   ArrheniusTe)
SELECT_RATE_LAWS(EXCHANGE,                   ArrheniusT,    ArrheniusT)
SELECT_RATE_LAWS(EXCITATION_M,               ArrheniusT,    ArrheniusT)
SELECT_RATE_LAWS(EXCITATION_E,               ArrheniusTe,   ArrheniusTe)

#undef SELECT_RATE_LAWS

//==============================================================================

RateManager::RateManager(size_t ns, const std::vector<Reaction>& reactions)
    : m_ns(ns), m_nr(reactions.size()), mp_lnkf(NULL), mp_lnkb(NULL),
      mp_gibbs(NULL)
{
    // Add all of the reactions' rate coefficients to the manager
    const size_t nr = reactions.size();
    for (size_t i = 0; i < m_nr; ++i)
        addReaction(i, reactions[i]);
    
    // Allocate storage in one block for both rate coefficient arrays and
    // species gibbs free energies 
    const size_t block_size = 2*m_nr + ns;
    mp_lnkf  = new double [block_size];
    mp_lnkb  = mp_lnkf + m_nr;
    mp_gibbs = mp_lnkb + m_nr;
    
    // Initialize the arrays to zero
    std::fill(mp_lnkf, mp_lnkf+block_size, 0.0);
}

//==============================================================================

RateManager::~RateManager()
{
    // Storage for lnkff and lnkfb were both allocated in mp_lnkff
    if (mp_lnkf != NULL)
        delete [] mp_lnkf;
}

//==============================================================================

void RateManager::addReaction(const size_t rxn, const Reaction& reaction)
{
    // Get the rate law which is being used in this reaction
    const RateLaw* p_rate = reaction.rateLaw();
    
    // Arrhenius reactions
    if (typeid(*p_rate) == typeid(Arrhenius)) {
        selectRate<MAX_REACTION_TYPES-1>(rxn, reaction);
    } else {
        throw InvalidInputError("rate law", typeid(*p_rate).name())
            << "Rate law is not implemented in RateManager.";
    }
}

//==============================================================================

template <int NReactionTypes>
void RateManager::selectRate(
    const size_t rxn, const Reaction& reaction)
{
    if (reaction.type() == NReactionTypes)
        addRate<
            typename RateSelector<NReactionTypes>::ForwardGroup,
            typename RateSelector<NReactionTypes>::ReverseGroup>(
            rxn, reaction);
    else
        selectRate<NReactionTypes-1>(rxn, reaction);
}

template <>
void RateManager::selectRate<0>(
    const size_t rxn, const Reaction& reaction)
{
    addRate<
        RateSelector<0>::ForwardGroup,
        RateSelector<0>::ReverseGroup>(rxn, reaction);
}

//==============================================================================

template <typename T, typename U>
struct is_same {
    enum { value = 0 };
};

template <typename T>
struct is_same<T,T> {
    enum { value = 1 };
};

//==============================================================================

template <typename ForwardGroup, typename ReverseGroup>
void RateManager::addRate(const size_t rxn, const Reaction& reaction)
{    
    m_rate_groups.addRateCoefficient<ForwardGroup>(rxn, reaction.rateLaw());
    
    if (reaction.isReversible()) {
        
        // Make use of forward computation when possible
        if (is_same<ForwardGroup, ReverseGroup>::value)
            m_to_copy.push_back(rxn);
        else
            // Evaluate at the reverse temperature
            // note: mp_lnkff+(rxn+m_nr) = mp_lnkfb+rxn
            m_rate_groups.addRateCoefficient<ReverseGroup>(
                rxn+m_nr, reaction.rateLaw());
        
        m_rate_groups.addReaction<ReverseGroup>(rxn, reaction);
        
    } else {
        m_irr.push_back(rxn);
    }
}

//==============================================================================

void RateManager::update(const Thermodynamics::Thermodynamics& thermo)
{
    // Evaluate all of the different rate coefficients
    m_rate_groups.logOfRateCoefficients(thermo.state(), mp_lnkf);
    
    // Copy rate coefficients which are the same as one of the previously
    // calculated ones
    std::vector<size_t>::const_iterator iter = m_to_copy.begin();
    for ( ; iter != m_to_copy.end(); ++iter) {
        const size_t index = *iter;
        mp_lnkb[index] = mp_lnkf[index];
    }
    
    // Subtract lnkeq(Tb) rate constants from the lnkf(Tb) to get lnkb(Tb)
    m_rate_groups.subtractLnKeq(thermo, mp_gibbs, mp_lnkb);
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation

