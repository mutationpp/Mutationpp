/**
 * @file StoichiometryManager.h
 *
 * Defines the Stoich<N> and StoichiometryManager classes which are heavily
 * based on the kinetics classes design from the Cantera kinetics library.
 */

#ifndef STOICHIOMETRYMANAGER_H
#define STOICHIOMETRYMANAGER_H

#include <vector>
#include "Numerics.h"


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
    
    /**
     * Performs the following operation on this reaction:
     * \f[
     *     r_j = r_j \prod_{i=1}^{n_s} s_i^{\nu_{ij}}
     * \f]
     * where \f$r_j\f$ is a reaction dependent quantity, \f$s_i\f$ are species
     * dependent quantities, and \f$\nu_{ij}\f$ is the stoichiometric coefficient
     * for species i in reaction j.
     */
    inline void multReaction(
        const Numerics::RealVector& s, Numerics::RealVector& r) const {
        for (int i = 0; i < N; ++i)
            r(m_rxn) *= s(m_sps[i]);
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
    inline void incrReaction(
        const Numerics::RealVector& s, Numerics::RealVector& r) const {
        for (int i = 0; i < N; ++i)
            r(m_rxn) += s(m_sps[i]);
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
    inline void decrReaction(
        const Numerics::RealVector& s, Numerics::RealVector& r) const {
        for (int i = 0; i < N; ++i)
            r(m_rxn) -= s(m_sps[i]);
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
    inline void incrSpecies(
        const Numerics::RealVector& r, Numerics::RealVector& s) const {
        for (int i = 0; i < N; ++i)
            s(m_sps[i]) += r(m_rxn);
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
    inline void decrSpecies(
        const Numerics::RealVector& r, Numerics::RealVector& s) const {
        for (int i = 0; i < N; ++i)
            s(m_sps[i]) -= r(m_rxn);
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
    
    void addReaction(const int rxn, const std::vector<size_t> &sps);
    
    void multReactions(
        const Numerics::RealVector& s, Numerics::RealVector& r) const;
    
    void incrReactions(
        const Numerics::RealVector& s, Numerics::RealVector& r) const;
    
    void decrReactions(
        const Numerics::RealVector& s, Numerics::RealVector& r) const;
    
    void incrSpecies(
        const Numerics::RealVector& r, Numerics::RealVector& s) const;
    
    void decrSpecies(
        const Numerics::RealVector& r, Numerics::RealVector& s) const;

private:

    std::vector<Stoich1> m_stoich1_vec;
    std::vector<Stoich2> m_stoich2_vec;
    std::vector<Stoich3> m_stoich3_vec;

}; // class StoichiometryManager

#endif // STOICHIOMETRYMANAGER_H
