#include "JacobianManager.h"
#include "Reaction.h"
#include "Numerics.h"

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
    
    for (int i = 0; i < Reactants::nSpecies(); ++i) {
        for (int j = 0; j < Reactants::nSpecies(); ++j)
            sjac[m_reacs(i)*ns + m_reacs(j)] -=
                m_reacs[i] * work[m_reacs(j)];
        
        for (int j = 0; j < Products::nSpecies(); ++j)
            sjac[m_reacs(i)*ns + m_prods(j)] -=
                m_reacs[i] * work[m_prods(j)];
    }
    
    for (int i = 0; i < Products::nSpecies(); ++i) {
        for (int j = 0; j < Reactants::nSpecies(); ++j)
            sjac[m_prods(i)*ns + m_reacs(j)] +=
                m_prods[i] * work[m_reacs(j)];
                
        for (int j = 0; j < Products::nSpecies(); ++j)
            sjac[m_prods(i)*ns + m_prods(j)] +=
                m_prods[i] * work[m_prods(j)];
    }
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
    
    /*cout << "kf: " << kf << endl;
    cout << "kb: " << kb << endl;
    cout << "RRf: " << rrf << endl;
    cout << "RRb: " << rrb << endl;
    cout << "RR: " << rr << endl;
    cout << "TB: " << tb << endl;*/

    m_reacs.diffRR(kf, conc, work, PlusEqualsTimes(tb));
    m_prods.diffRR(kb, conc, work, MinusEqualsTimes(tb));
    
    //cout << "reactants: ";
    for (int i = 0; i < Reactants::nSpecies(); ++i) {
        //cout << setw(3) << m_reacs(i);
        for (int j = 0; j < ns; ++j)
            sjac[m_reacs(i)*ns + j] -= m_reacs[i] * work[j];
    }
    
    //cout << endl << "products: ";
    for (int i = 0; i < Products::nSpecies(); ++i) {
        //cout << setw(3) << m_prods(i);
        for (int j = 0; j < ns; ++j)
            sjac[m_prods(i)*ns + j] += m_prods[i] * work[j];
    }
    //cout << endl;
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
        std::cout << "Error: reactants' stoichiometry is invalid for reaction"
             << "\"" << reaction.formula() << "\"!" << std::endl;
        exit(1);
    }
    
    if (!getJacStoich(reaction.products(), &p_prods, prodsType)) {
        std::cout << "Error: products' stoichiometry is invalid for reaction"
             << "\"" << reaction.formula() << "\"!" << std::endl;
        exit(1);
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

void JacobianManager::computeJacobian(
    const double* const kf, const double* const kb, const double* const conc, 
    double* const sjac) const
{
    const size_t ns = m_thermo.nSpecies();
    const size_t nr = m_reactions.size();
    
    // Make sure we are staring with a clean slate
    for (int i = 0; i < ns*ns; ++i)
        sjac[i] = 0.0;

    // Loop over each reaction and compute the dRR/dconc_k
    for (int i = 0; i < nr; ++i)
        m_reactions[i]->contributeToJacobian(
            kf[i], kb[i], conc, mp_work, sjac, ns);
    
    // Finally, multiply by the species molecular weight ratios
    for (int i = 0; i < ns; ++i)
        mp_work[i] = m_thermo.speciesMw(i);
    
    for (int i = 0, index = 0; i < ns; ++i) {
        for (int j = 0; j < ns; ++j, ++index)
            sjac[index] *= mp_work[i] / mp_work[j];
    }
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation

