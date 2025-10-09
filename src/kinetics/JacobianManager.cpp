/**
 * @file JacobianManager.cpp
 *
 * @brief Implements Kinetics::JacobianManager class and other helper classes.
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


#include "JacobianManager.h"
#include "Reaction.h"
#include "Functors.h"

#include <iomanip>
#include <iostream>

typedef Mutation::Numerics::Equals<double>            Equals;
typedef Mutation::Numerics::MinusEquals<double>       MinusEquals;
typedef Mutation::Numerics::PlusEqualsTimes<double>   PlusEqualsTimes;
typedef Mutation::Numerics::MinusEqualsTimes<double>  MinusEqualsTimes;

namespace Mutation {
    namespace Kinetics {

//==============================================================================

template <typename Reactants, typename Products>
void swap(
    ThirdbodyReactionStoich<Reactants, Products>& left,
    ThirdbodyReactionStoich<Reactants, Products>& right)
{
    std::swap(left.m_reacs,  right.m_reacs);
    std::swap(left.m_prods,  right.m_prods);
    std::swap(left.mp_alpha, right.mp_alpha);
    std::swap(left.m_ns,     right.m_ns);
}

//==============================================================================

template <typename Reactants, typename Products>
void ReactionStoich<Reactants, Products>::contributeToJacobian(
    const double kf, const double kb, const double* const conc, 
    double* const work, double* const sjac, const size_t ns) const
{
    // Need to make sure that product species are zeroed out in the work
    // array (don't need to zero out whole array)
    for (int i = 0; i < Products::nSpecies(); ++i)
        work[m_prods(i)] = 0.0;

    // Compute the derivative of reaction rate with respect to the reactants
    // and products (all other terms are zero)
    m_reacs.diffRR(kf, conc, work, Equals());
    m_prods.diffRR(kb, conc, work, MinusEquals());

    // Only loop over the necessary species
    for (auto& pi : m_index_stoich)
        for (auto& pj : m_index_stoich)
            sjac[pi.first*ns+pj.first] += pi.second * work[pj.first];
}

//==============================================================================

template <typename Reactants, typename Products>
void ReactionStoich<Reactants, Products>::contributeToNetProgressJacobian(
    const double kf, const double kb, const double* const conc, 
    double* const work, double* const p_dropRho, const size_t ns, const size_t r) const
{
    // Need to make sure that product species are zeroed out in the work
    // array (don't need to zero out whole array)
    for (int i = 0; i < ns; ++i)
        work[i] = 0.0;

    // Compute the derivative of reaction rate with respect to the reactants
    // and products (all other terms are zero)
    m_reacs.diffRR(kf, conc, work, Equals());
    m_prods.diffRR(kb, conc, work, MinusEquals());

    // Only loop over the necessary species
    for (int i = 0; i < ns; ++i)
            p_dropRho[i+r*ns] = work[i];

}

//==============================================================================

template <typename Reactants, typename Products>
void ThirdbodyReactionStoich<Reactants, Products>::contributeToJacobian(
    const double kf, const double kb, const double* const conc, 
    double* const work, double* const sjac, const size_t ns) const
{
    const double rrf = m_reacs.rr(kf, conc);
    const double rrb = m_prods.rr(kb, conc);
    const double rr = rrf - rrb;
    double tb = 0.0;
    
    for (int i = 0; i < ns; ++i) {
        work[i] = mp_alpha[i] * rr;
        tb += mp_alpha[i] * conc[i];
    }

    m_reacs.diffRR(kf, conc, work, PlusEqualsTimes(tb));
    m_prods.diffRR(kb, conc, work, MinusEqualsTimes(tb));

    for (auto& p : m_index_stoich)
        for (int j = 0; j < ns; ++j)
            sjac[p.first*ns+j] += p.second * work[j];
}

//==============================================================================

template <typename Reactants, typename Products>
void ThirdbodyReactionStoich<Reactants, Products>::contributeToNetProgressJacobian(
    const double kf, const double kb, const double* const conc, 
    double* const work, double* const p_dropRho, const size_t ns, const size_t r) const
{
    const double rrf = m_reacs.rr(kf, conc);
    const double rrb = m_prods.rr(kb, conc);
    const double rr = rrf - rrb;
    double tb = 0.0;
    
    for (int i = 0; i < ns; ++i)
        work[i] = 0.0;

    for (int i = 0; i < ns; ++i) {
        work[i] = mp_alpha[i] * rr;
        tb += mp_alpha[i] * conc[i];
    }

    m_reacs.diffRR(kf, conc, work, PlusEqualsTimes(tb));
    m_prods.diffRR(kb, conc, work, MinusEqualsTimes(tb));

    for (int i = 0; i < ns; ++i)
            p_dropRho[i+r*ns] = work[i];

}

//==============================================================================

template <typename Reactants>
void JacobianManager::addReactionStoich(
    JacStoichBase* p_reacs, JacStoichBase* p_prods, const StoichType type, 
    const Reaction& reaction)
{
    if (reaction.isThirdbody()) {
        // Thirdbody reactions        
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_work[i] = 1.0;
        if (m_thermo.hasElectrons())
            mp_work[0] = 0.0;
        
        for (int i = 0; i < reaction.efficiencies().size(); ++i)
            mp_work[reaction.efficiencies()[i].first] = 
                reaction.efficiencies()[i].second;
        
        switch (type) {
            case STOICH_11:
                m_reactions.push_back(
                    new ThirdbodyReactionStoich<Reactants, JacStoich11>(
                        p_reacs, p_prods, mp_work, m_thermo.nSpecies()));
                break;
            case STOICH_21:
                m_reactions.push_back(
                    new ThirdbodyReactionStoich<Reactants, JacStoich21>(
                        p_reacs, p_prods, mp_work, m_thermo.nSpecies()));
                break;
            case STOICH_22:
                m_reactions.push_back(
                    new ThirdbodyReactionStoich<Reactants, JacStoich22>(
                        p_reacs, p_prods, mp_work, m_thermo.nSpecies()));
                break;
            case STOICH_31:
                m_reactions.push_back(
                    new ThirdbodyReactionStoich<Reactants, JacStoich31>(
                        p_reacs, p_prods, mp_work, m_thermo.nSpecies()));
                break;
            case STOICH_32:
                m_reactions.push_back(
                    new ThirdbodyReactionStoich<Reactants, JacStoich32>(
                        p_reacs, p_prods, mp_work, m_thermo.nSpecies()));
                break;
            case STOICH_33:
                m_reactions.push_back(
                    new ThirdbodyReactionStoich<Reactants, JacStoich33>(
                        p_reacs, p_prods, mp_work, m_thermo.nSpecies()));
                break;
        }
    } else {
        // Not thirdbody reactions
        switch (type) {
            case STOICH_11:
                m_reactions.push_back(
                    new ReactionStoich<Reactants, JacStoich11>(
                        p_reacs, p_prods));
                break;
            case STOICH_21:
                m_reactions.push_back(
                    new ReactionStoich<Reactants, JacStoich21>(
                        p_reacs, p_prods));
                break;
            case STOICH_22:
                m_reactions.push_back(
                    new ReactionStoich<Reactants, JacStoich22>(
                        p_reacs, p_prods));
                break;
            case STOICH_31:
                m_reactions.push_back(
                    new ReactionStoich<Reactants, JacStoich31>(
                        p_reacs, p_prods));
                break;
            case STOICH_32:
                m_reactions.push_back(
                    new ReactionStoich<Reactants, JacStoich32>(
                        p_reacs, p_prods));
                break;
            case STOICH_33:
                m_reactions.push_back(
                    new ReactionStoich<Reactants, JacStoich33>(
                        p_reacs, p_prods));
                break;
        }
    }
}

//==============================================================================

bool JacobianManager::getJacStoich(
    const std::vector<int>& stoich_vec, JacStoichBase** p_stoich,
    StoichType& type) const
{    
    switch (stoich_vec.size()) {
        case 1:
            *p_stoich = new JacStoich11(stoich_vec[0]);
            type = STOICH_11;
            return true;
        case 2:
            if (stoich_vec[0] == stoich_vec[1]) {
               *p_stoich = new JacStoich21(stoich_vec[0]);
                type = STOICH_21;
            } else {
                *p_stoich = new JacStoich22(stoich_vec[0], stoich_vec[1]);
                type = STOICH_22;
            }
            return true;
        case 3:
            if (stoich_vec[0] == stoich_vec[1]) {
                if (stoich_vec[1] == stoich_vec[2]) {
                    *p_stoich = new JacStoich31(stoich_vec[0]);
                    type = STOICH_31;
                } else {
                    *p_stoich = new JacStoich32(stoich_vec[0], stoich_vec[2]);
                    type = STOICH_32;
                }
            } else {
                if (stoich_vec[1] == stoich_vec[2]) {
                    *p_stoich = new JacStoich32(stoich_vec[1], stoich_vec[0]);
                    type = STOICH_32;
                } else if (stoich_vec[0] == stoich_vec[2]) {
                    *p_stoich = new JacStoich32(stoich_vec[0], stoich_vec[1]);
                    type = STOICH_32;
                } else {
                    *p_stoich = new JacStoich33(
                        stoich_vec[0], stoich_vec[1], stoich_vec[2]);
                    type = STOICH_33;
                }
            }
            return true;
        default:
            return false;
    }
}

//==============================================================================

void JacobianManager::addReaction(const Reaction& reaction)
{
    // Need to figure out which class of reaction stoichiometries this 
    // reaction belongs to
    JacStoichBase* p_reacs;
    JacStoichBase* p_prods;
    StoichType reacsType;
    StoichType prodsType;
    
    if (!getJacStoich(reaction.reactants(), &p_reacs, reacsType)) {
        throw InvalidInputError("reaction", reaction.formula())
            << "Reactants' stoichiometry is not implemented in "
            << "JacobianManager.";
    }
    
    if (!getJacStoich(reaction.products(), &p_prods, prodsType)) {
        throw InvalidInputError("reaction", reaction.formula())
            << "Products' stoichiometry is not implemented in "
            << "JacobianManager.";
    }
    
    switch (reacsType) {
        case STOICH_11:
            addReactionStoich<JacStoich11>(
                p_reacs, p_prods, prodsType, reaction);
            break;
        case STOICH_21:
            addReactionStoich<JacStoich21>(
                p_reacs, p_prods, prodsType, reaction);
            break;
        case STOICH_22:
            addReactionStoich<JacStoich22>(
                p_reacs, p_prods, prodsType, reaction);
            break;
        case STOICH_31:
            addReactionStoich<JacStoich31>(
                p_reacs, p_prods, prodsType, reaction);
            break;
        case STOICH_32:
            addReactionStoich<JacStoich32>(
                p_reacs, p_prods, prodsType, reaction);
            break;
        case STOICH_33:
            addReactionStoich<JacStoich33>(
                p_reacs, p_prods, prodsType, reaction);
            break;
    }
    
    delete p_reacs;
    delete p_prods;
}

//==============================================================================

void JacobianManager::computeTReactionTDerivatives(
     const std::vector<Reaction>& reactions, const double* const p_t,
     std::vector<std::pair<double, double> >& dTrxndT) const
    {
        
        for(size_t i = 0; i < reactions.size(); ++i) {
            const size_t type = reactions[i].type();
            switch (type) {
                case ASSOCIATIVE_IONIZATION:
                    dTrxndT.push_back({1.0, 0.0});
                    break;
                case DISSOCIATIVE_RECOMBINATION:
                    dTrxndT.push_back({0.0, 1.0});
                    break;
                case ASSOCIATIVE_DETACHMENT:
                    dTrxndT.push_back({1.0, 0.0});
                    break;
                case DISSOCIATIVE_ATTACHMENT:
                    dTrxndT.push_back({0.0, 1.0});
                    break;
                case DISSOCIATION_M:
                    dTrxndT.push_back({0.5*sqrt(p_t[0]*p_t[1])/p_t[0], 1.0});
                    break;
                case DISSOCIATION_E:
                    dTrxndT.push_back({0.0, 0.0});
                    break;
                case RECOMBINATION_M:
                    dTrxndT.push_back({1.0, 0.5*sqrt(p_t[0]*p_t[1])/p_t[0]});
                    break;
                case RECOMBINATION_E:
                    dTrxndT.push_back({0.0, 0.0});
                    break;
                case IONIZATION_M:
                    dTrxndT.push_back({1.0, 1.0});
                    break;
                case IONIZATION_E:
                    dTrxndT.push_back({0.0, 1.0});
                    break;
                case ION_RECOMBINATION_M:
                    dTrxndT.push_back({1.0, 1.0});
                    break;
                case ION_RECOMBINATION_E:
                    dTrxndT.push_back({1.0, 0.0});
                    break;
                case ELECTRONIC_ATTACHMENT_E:
                    dTrxndT.push_back({0.0, 0.0});
                    break;
                case ELECTRONIC_DETACHMENT_E:
                    dTrxndT.push_back({0.0, 0.0});
                    break;
                case ELECTRONIC_ATTACHMENT_M:
                    dTrxndT.push_back({0.0, 1.0});
                    break;
                case ELECTRONIC_DETACHMENT_M:
                    dTrxndT.push_back({1.0, 0.0});
                    break;
                case EXCHANGE:
                    dTrxndT.push_back({1.0, 1.0});
                    break;
                case EXCITATION_E:
                    dTrxndT.push_back({0.0, 0.0});
                    break;
                case EXCITATION_M:
                    dTrxndT.push_back({1.0, 1.0});
                    break;
                default:
                    dTrxndT.push_back({0.0, 0.0});
                    break;
            }
        }
    }

//==============================================================================

void JacobianManager::computeTReactionTvDerivatives(
     const std::vector<Reaction>& reactions, const double* const p_t,
     std::vector<std::pair<double, double> >& dTrxndTv) const
    {
        
        for(size_t i = 0; i < reactions.size(); ++i) {
            const size_t type = reactions[i].type();
            switch (type) {
                case ASSOCIATIVE_IONIZATION:
                    dTrxndTv.push_back({0.0, 1.0});
                    break;
                case DISSOCIATIVE_RECOMBINATION:
                    dTrxndTv.push_back({1.0, 0.0});
                    break;
                case ASSOCIATIVE_DETACHMENT:
                    dTrxndTv.push_back({0.0, 1.0});
                    break;
                case DISSOCIATIVE_ATTACHMENT:
                    dTrxndTv.push_back({1.0, 0.0});
                    break;
                case DISSOCIATION_M:
                    dTrxndTv.push_back({0.5*sqrt(p_t[0]*p_t[1])/p_t[1], 0.0});
                    break;
                case DISSOCIATION_E:
                    dTrxndTv.push_back({1.0, 1.0});
                    break;
                case RECOMBINATION_M:
                    dTrxndTv.push_back({0.0, 0.5*sqrt(p_t[0]*p_t[1])/p_t[1]});
                    break;
                case RECOMBINATION_E:
                    dTrxndTv.push_back({1.0, 1.0});
                    break;
                case IONIZATION_M:
                    dTrxndTv.push_back({0.0, 0.0});
                    break;
                case IONIZATION_E:
                    dTrxndTv.push_back({1.0, 0.0});
                    break;
                case ION_RECOMBINATION_M:
                    dTrxndTv.push_back({0.0, 0.0});
                    break;
                case ION_RECOMBINATION_E:
                    dTrxndTv.push_back({0.0, 1.0});
                    break;
                case ELECTRONIC_ATTACHMENT_E:
                    dTrxndTv.push_back({1.0, 1.0});
                    break;
                case ELECTRONIC_DETACHMENT_E:
                    dTrxndTv.push_back({1.0, 1.0});
                    break;
                case ELECTRONIC_ATTACHMENT_M:
                    dTrxndTv.push_back({1.0, 0.0});
                    break;
                case ELECTRONIC_DETACHMENT_M:
                    dTrxndTv.push_back({0.0, 1.0});
                    break;
                case EXCHANGE:
                    dTrxndTv.push_back({0.0, 0.0});
                    break;
                case EXCITATION_E:
                    dTrxndTv.push_back({1.0, 1.0});
                    break;
                case EXCITATION_M:
                    dTrxndTv.push_back({0.0, 0.0});
                    break;
                default:
                    dTrxndTv.push_back({0.0, 0.0});
                    break;
            }
        }
    }

//==============================================================================

void JacobianManager::computeNetRateProgressJacobian(
    const double* const kf, const double* const kb, const double* const conc, 
    double* const p_dropRho) const
{
    const size_t ns = m_thermo.nSpecies();
    const size_t nr = m_reactions.size();
    
    // Make sure we are staring with a clean slate
    std::fill(p_dropRho, p_dropRho + ns*nr, 0.0);

    // Loop over each reaction and compute the dRR/dconc_k
    for (int i = 0; i < nr; ++i) 
        m_reactions[i]->contributeToNetProgressJacobian(
            kf[i], kb[i], conc, mp_work1, p_dropRho, ns, i);
    
    
    // Finally, multiply by the species molecular weight ratios
    for (int i = 0, index = 0; i < nr; ++i) 
        for (int j = 0; j < ns; ++j, ++index)
            p_dropRho[index] /= m_thermo.speciesMw(j);
    
}

//==============================================================================

void JacobianManager::computeJacobian(
    const double* const kf, const double* const kb, const double* const conc, 
    double* const sjac) const
{
    const size_t ns = m_thermo.nSpecies();
    const size_t nr = m_reactions.size();
    
    // Make sure we are staring with a clean slate
    std::fill(sjac, sjac+ns*ns, 0.0);

    // Loop over each reaction and compute the dRR/dconc_k
    for (int i = 0; i < nr; ++i) 
        m_reactions[i]->contributeToJacobian(
            kf[i], kb[i], conc, mp_work, sjac, ns);
    
    
    // Finally, multiply by the species molecular weight ratios
    for (int i = 0, index = 0; i < ns; ++i) 
        for (int j = 0; j < ns; ++j, ++index)
            sjac[index] *= m_thermo.speciesMw(i) / m_thermo.speciesMw(j);
}

//==============================================================================

void JacobianManager::computeJacobianT(
    const double* const p_dropf, const double* const p_dropb, 
    const std::vector<Reaction>& reactions, double* const p_dropT) const
    {
        const size_t nt = m_thermo.nEnergyEqns();
        std::vector<std::pair<double, double> > dTrxndT;
        
        std::fill(mp_T, mp_T + nt, 0.);
        m_thermo.getTemperatures(mp_T);
        
        computeTReactionTDerivatives(reactions, mp_T, dTrxndT);
        
        for (int i = 0; i < reactions.size(); ++i)
            p_dropT[i] = p_dropf[i]*dTrxndT[i].first - p_dropb[i]*dTrxndT[i].second;
    }

//==============================================================================

void JacobianManager::computeJacobianTv(
    const double* const p_dropf, const double* const p_dropb, 
    const std::vector<Reaction>& reactions, double* const p_dropTv) const
    {
        const size_t nt = m_thermo.nEnergyEqns();
        std::vector<std::pair<double, double> > dTrxndTv;
        
        std::fill(mp_T, mp_T + nt, 0.);
        m_thermo.getTemperatures(mp_T);
        
        computeTReactionTvDerivatives(reactions, mp_T, dTrxndTv);
        
        for (int i = 0; i < reactions.size(); ++i)
            p_dropTv[i] = p_dropf[i]*dTrxndTv[i].first - p_dropb[i]*dTrxndTv[i].second;
    }

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation

