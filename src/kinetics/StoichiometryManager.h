/**
 * @file StoichiometryManager.h
 *
 * @brief Defines the Stoich<N> and StoichiometryManager classes which are
 * inspired from the kinetics classes of the Cantera library.
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

#ifndef KINETICS_STOICHIOMETRY_MANAGER_H
#define KINETICS_STOICHIOMETRY_MANAGER_H

#include <vector>
#include <cstdlib>

namespace Mutation {
    namespace Kinetics {


/**
 * Represents a single reaction's stoichiometry.
 *
 * In fact, this class represents one side of a complete reaction.  Its purpose
 * is to provide a sparse-matrix type functionality for computing stoichiometry
 * dependent quantities by storing the dependencies of reactions on the species
 * involved in them.  This is done in the form of species and reaction indices
 * used in global arrays. This class is templated based on the number of species
 * present in the stoichiometry description.
 */
template <int N>
class Stoich
{
public:

    /**
     * Constructor.
     */
    Stoich()
    {
        m_rxn = 0;
        for (int i = 0; i < N; ++i)
            m_sps[i] = 0;
    }
    
    /*template <typename OP>
    inline void updateReaction(
        const double* const p_s, double* const p_r, const OP& op) const
    {
        for (int i = 0; i < N; ++i)
            op(p_r[m_rnx], p_s[m_sps[i]]);
    }*/
    
    /**
     * Performs the following operation on this reaction:
     * \f[
     *     r_j = r_j \prod_{i=1}^{n_s} s_i^{\nu_{ij}}
     * \f]
     * where \f$r_j\f$ is a reaction dependent quantity, \f$s_i\f$ are species
     * dependent quantities, and \f$\nu_{ij}\f$ is the stoichiometric coefficient
     * for species i in reaction j.
     */
    inline void multReaction(const double* const p_s, double* const p_r) const
    {
        for (int i = 0; i < N; ++i)
            p_r[m_rxn] *= p_s[m_sps[i]];
    }
    
    /**
     * Performs the following operation on this reaction:
     * \f[
     *     r_j = r_j + \sum_{i=1}^{n_s} \nu_{ij} s_i
     * \f]
     * where \f$r_j\f$ is a reaction dependent quantity, \f$s_i\f$ are species
     * dependent quantities, and \f$\nu_{ij}\f$ is the stoichiometric coefficient
     * for species i in reaction j.
     */
    inline void incrReaction(const double* const p_s, double* const p_r) const
    {
        for (int i = 0; i < N; ++i)
            p_r[m_rxn] += p_s[m_sps[i]];
    }
    
    /**
     * Performs the following operation on this reaction:
     * \f[
     *     r_j = r_j - \sum_{i=1}^{n_s} \nu_{ij} s_i
     * \f]
     * where \f$r_j\f$ is a reaction dependent quantity, \f$s_i\f$ are species
     * dependent quantities, and \f$\nu_{ij}\f$ is the stoichiometric coefficient
     * for species i in reaction j.
     */
    inline void decrReaction(const double* const p_s, double* const p_r) const 
    {
        for (int i = 0; i < N; ++i)
            p_r[m_rxn] -= p_s[m_sps[i]];
    }
    
    /**
     * Performs the following operation on the species in this reaction:
     * \f[
     *     s_i = s_i + \nu_{ij} r_j, \quad \forall \; i = 1,\dots,n_s
     * \f]
     * where \f$r_j\f$ is a reaction dependent quantity, \f$s_i\f$ are species
     * dependent quantities, and \f$\nu_{ij}\f$ is the stoichiometric coefficient
     * for species i in reaction j.
     */
    inline void incrSpecies(const double* const p_r, double* const p_s) const
    {
        for (int i = 0; i < N; ++i)
            p_s[m_sps[i]] += p_r[m_rxn];
    }
    
    /**
     * Performs the following operation on the species in this reaction:
     * \f[
     *     s_i = s_i - \nu_{ij} r_j, \quad \forall \; i = 1,\dots,n_s
     * \f]
     * where \f$r_j\f$ is a reaction dependent quantity, \f$s_i\f$ are species
     * dependent quantities, and \f$\nu_{ij}\f$ is the stoichiometric coefficient
     * for species i in reaction j.
     */
    inline void decrSpecies(const double* const p_r, double* const p_s) const
    {
        for (int i = 0; i < N; ++i)
            p_s[m_sps[i]] -= p_r[m_rxn];
    }
    
protected:

    size_t m_rxn;
    size_t m_sps[N];
    
}; // class Stoich<N>


/**
 * Extends the general Stoich base class for a single reactant.
 */
class Stoich1 : public Stoich<1>
{
public:
    Stoich1(size_t rxn, size_t sp1)
    { 
        m_rxn = rxn;
        m_sps[0] = sp1;
    }
};

/**
 * Extends the general Stoich base class for two reactants.
 */
class Stoich2 : public Stoich<2>
{
public:
    Stoich2(size_t rxn, size_t sp1, size_t sp2)
    { 
        m_rxn = rxn;
        if (sp1 <= sp2) {
            m_sps[0] = sp1;
            m_sps[1] = sp2;
        } else {
            m_sps[0] = sp2;
            m_sps[1] = sp1;
        }
    }
};


/**
 * Extends the general Stoich base class for three reactants.
 */
class Stoich3 : public Stoich<3>
{
public:
    Stoich3(size_t rxn, size_t sp1, size_t sp2, size_t sp3)
    { 
        m_rxn = rxn;
        
        if (sp1 <= sp2) {
            m_sps[0] = sp1;
            m_sps[1] = sp2;
        } else {
            m_sps[0] = sp2;
            m_sps[1] = sp1;
        }
        
        if (m_sps[1] <= sp3)
            m_sps[2] = sp3;
        else {
            if (sp3 <= m_sps[0]) {
                m_sps[2] = m_sps[1];
                m_sps[1] = m_sps[0];
                m_sps[0] = sp3;
            } else {
                m_sps[2] = m_sps[1];
                m_sps[1] = sp3;
            }
        }
    }
};


/**
 * Provides sparse-matrix type operations for quantities that depend on reaction
 * stoichiometries.
 */
class StoichiometryManager
{
public:

    StoichiometryManager() { }
    
    void addReaction(const int rxn, const std::vector<int>& sps);
    
    void multReactions(const double* const p_s, double* const p_r) const;
    
    void incrReactions(const double* const p_s, double* const p_r) const;
    
    void decrReactions(const double* const p_s, double* const p_r) const;
    
    void incrSpecies(const double* const p_r, double* const p_s) const;
    
    void decrSpecies(const double* const p_r, double* const p_s) const;

private:

    std::vector<Stoich1> m_stoich1_vec;
    std::vector<Stoich2> m_stoich2_vec;
    std::vector<Stoich3> m_stoich3_vec;

}; // class StoichiometryManager


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_STOICHIOMETRY_MANAGER_H
