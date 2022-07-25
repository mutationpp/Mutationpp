/**
 * @file JacobianManager.h
 *
 * @brief Declaration of JacobianManagar class and its helper classes.
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


#ifndef JACOBIAN_MANAGER_H
#define JACOBIAN_MANAGER_H

#include "Thermodynamics.h"
#include <iostream>

namespace Mutation {
    namespace Kinetics {

class Reaction;

/**
 * Enumerates the different JacStoich** types.
 */
enum StoichType {
    STOICH_11,
    STOICH_21,
    STOICH_22,
    STOICH_31,
    STOICH_32,
    STOICH_33
};

/**
 * Empty base class for all JacStoich** types in order to facilitate 
 * initialization of the JacobianManager class.
 */
class JacStoichBase
{
public:
    virtual ~JacStoichBase() {};
};

/**
 * Represents "A" stoichiometry.
 */
class JacStoich11 : public JacStoichBase
{
public:
    
    JacStoich11(const size_t sp)
        : m_sp(sp)
    { }

    static int nSpecies() {
        return 1;
    }
    
    int operator () (const int& index) const {
        return m_sp;
    }
    
    double operator [] (const int& index) const {
        return 1.0;
    }
    
    inline double rr(double k, const double* const conc) const {
        return k * conc[m_sp];
    }

    template <typename Op>
    inline void diffRR(
        const double& k, const double* const conc, double* const drr, 
        const Op& op) const 
    {
        op(drr[m_sp], k);
    }

private:

    size_t m_sp;

};

/**
 * Represents "2A" stoichiometry.
 */
class JacStoich21 : public JacStoichBase
{
public:
    
    JacStoich21(const size_t sp)
        : m_sp(sp)
    { }

    static int nSpecies() {
        return 1;
    }
    
    int operator () (const int& index) const {
        return m_sp;
    }
    
    double operator [] (const int& index) const {
        return 2.0;
    }

    inline double rr(double k, const double* const conc) const {
        return k * conc[m_sp] * conc[m_sp];
    }

    template <typename Op>
    inline void diffRR(
        const double& k, const double* const conc, double* const drr, 
        const Op& op) const 
    {
        op(drr[m_sp], 2.0 * k * conc[m_sp]);
    }

private:
    
    size_t m_sp;
    
};

/**
 * Represents "A+B" stoichiometry.
 */
class JacStoich22 : public JacStoichBase
{
public:
    
    JacStoich22(const size_t sp1, const size_t sp2)
        : m_sp1(sp1), m_sp2(sp2)
    { }

    static int nSpecies() {
        return 2;
    }
    
    int operator () (const int& index) const {
        return (index == 0 ? m_sp1 : m_sp2);
    }
    
    double operator [] (const int& index) const {
        return 1.0;
    }
    
    inline double rr(double k, const double* const conc) const {
        return k * conc[m_sp1] * conc[m_sp2];
    }

    template <typename Op>
    inline void diffRR(
        const double k, const double* const conc, double* const drr, 
        const Op& op) const 
    {
        op(drr[m_sp1], k * conc[m_sp2]);
        op(drr[m_sp2], k * conc[m_sp1]);
    }

private:

    size_t m_sp1;
    size_t m_sp2;

};

/**
 * Represents "3A" stoichiometry.
 */
class JacStoich31 : public JacStoichBase
{
public:
    
    JacStoich31(const size_t sp)
        : m_sp(sp)
    { }

    static int nSpecies() {
        return 1;
    }
    
    int operator () (const int& index) const {
        return m_sp;
    }
    
    double operator [] (const int& index) const {
        return 3.0;
    }
    
    inline double rr(double k, const double* const conc) const {
        return k * conc[m_sp] * conc[m_sp] * conc[m_sp];
    }

    template <typename Op>
    inline void diffRR(
        const double k, const double* const conc, double* const drr, 
        const Op& op) const 
    {
        op(drr[m_sp], 3.0 * k * conc[m_sp] * conc[m_sp]);
    }
private:

    size_t m_sp;

};

/**
 * Represents "2A+B" stoichiometry.
 */
class JacStoich32 : public JacStoichBase
{
public:
    
    JacStoich32(const size_t sp1, const size_t sp2)
        : m_sp1(sp1), m_sp2(sp2)
    { }

    static int nSpecies() {
        return 2;
    }
    
    int operator () (const int& index) const {
        return (index == 0 ? m_sp1 : m_sp2);
    }
    
    double operator [] (const int& index) const {
        return (index == 0 ? 2.0 : 1.0);
    }
    
    inline double rr(double k, const double* const conc) const {
        return k * conc[m_sp1] * conc[m_sp1] * conc[m_sp2];
    }

    template <typename Op>
    inline void diffRR(
        const double k, const double* const conc, double* const drr, 
        const Op& op) const 
    {
        op(drr[m_sp1], 2.0 * k * conc[m_sp1] * conc[m_sp2]);
        op(drr[m_sp2], k * conc[m_sp1] * conc[m_sp1]);
    }

private:

    size_t m_sp1;
    size_t m_sp2;

};

/**
 * Represents "A+B+C" stoichiometry.
 */
class JacStoich33 : public JacStoichBase
{
public:
    
    JacStoich33(const size_t sp1, const size_t sp2, const size_t sp3)
        : m_sp1(sp1), m_sp2(sp2), m_sp3(sp3)
    { }

    static int nSpecies() {
        return 3;
    }
    
    int operator () (const int& index) const {
        return (index == 0 ? m_sp1 : (index == 1 ? m_sp2 : m_sp3));
    }
    
    double operator [] (const int& index) const {
        return 1.0;
    }
    
    inline double rr(double k, const double* const conc) const {
        return k * conc[m_sp1] * conc[m_sp2] * conc[m_sp3];
    }

    template <typename Op>
    inline void diffRR(
        const double k, const double* const conc, double* const drr, 
        const Op& op) const 
    {
        op(drr[m_sp1], k * conc[m_sp2] * conc[m_sp3]);
        op(drr[m_sp2], k * conc[m_sp1] * conc[m_sp3]);
        op(drr[m_sp3], k * conc[m_sp1] * conc[m_sp2]);
    }

private:

    size_t m_sp1;
    size_t m_sp2;
    size_t m_sp3;

};


/**
 * Abstract base class for all derived ReactionStoich types to enable
 * polymorphism.
 */
class ReactionStoichBase
{
public:
    
    /**
     * Destructor.
     */
    virtual ~ReactionStoichBase() { }

    virtual void contributeToJacobian(
        const double kf, const double kb, const double* const conc, 
        double* const work, double* const sjac, const size_t ns) const = 0;
};

/**
 * Connects reactant and product stoichiometry types together so that Jacobian
 * contributions made by a single reaction can be determined.
 */
template <typename Reactants, typename Products>
class ReactionStoich : public ReactionStoichBase
{
public:
    
    /**
     * Constructor.
     */
    ReactionStoich(JacStoichBase* reacs, JacStoichBase* prods)
        : m_reacs(*static_cast<Reactants*>(reacs)), 
          m_prods(*static_cast<Products*>(prods))
    { 
        // Save indices and reaction stoichiometry of the involved species
        for (int i = 0; i < Reactants::nSpecies(); ++i)
            m_index_stoich.emplace_back(m_reacs(i), -m_reacs[i]);
        
        for (int i = 0; i < Products::nSpecies(); ++i) {
            int k = -1;
            for (int j = 0; j < Reactants::nSpecies(); ++j)
                if (m_index_stoich[j].first == m_prods(i)) k = j;
            
            if (k == -1)
                m_index_stoich.emplace_back(m_prods(i), m_prods[i]);
            else
                m_index_stoich[k].second += m_prods[i];
        }
    }
    
    /**
     * Destructor.
     */
    virtual ~ReactionStoich() { }
    
    /**
     * Adds this reaction's contribution to the species Jacobian matrix.  Note
     * that this will affect the rows of each species that particpates in this
     * reaction.  The compiler should be able to unroll each loop since the
     * sizes are known at compile time given the Reactants and Products 
     * stoichiometry types.
     */
    void contributeToJacobian(
        const double kf, const double kb, const double* const conc, 
        double* const work, double* const sjac, const size_t ns) const;

protected:
    
    Reactants m_reacs;
    Products  m_prods;
    std::vector< std::pair<int, int> > m_index_stoich;
};


// Forward declaration of class ThirdbodyReactionStoich for swap.
template <typename Reactants, typename Products>
class ThirdbodyReactionStoich;

/**
 * Allows swapping of ThirdbodyReactionStoich.
 */    
template <typename Reactants, typename Products>
void swap(
    ThirdbodyReactionStoich<Reactants, Products>& left,
    ThirdbodyReactionStoich<Reactants, Products>& right);


/**
 * Connects reactant and product stoichiometry types together for a thirdbody
 * reaction so that Jacobian contributions made by a single reaction can be 
 * determined.
 */
template <typename Reactants, typename Products>
class ThirdbodyReactionStoich : public ReactionStoich<Reactants, Products>
{
public:

    using ReactionStoich<Reactants, Products>::m_reacs;
    using ReactionStoich<Reactants, Products>::m_prods;
    using ReactionStoich<Reactants, Products>::m_index_stoich;

    /**
     * Constructor.
     */
    ThirdbodyReactionStoich(
        JacStoichBase* reacs, JacStoichBase* prods, const double* const alpha,
        const size_t ns)
        : ReactionStoich<Reactants, Products>(reacs, prods), m_ns(ns),
          mp_alpha(new double [ns])
    { 
        std::copy(alpha, alpha+ns, mp_alpha);
    }
    
    /**
     * Copy constructor.
     */
    ThirdbodyReactionStoich(const ThirdbodyReactionStoich& to_copy)
        : ReactionStoich<Reactants, Products>(m_reacs, m_prods), 
          m_ns(to_copy.m_ns),
          mp_alpha(new double [to_copy.m_ns])
    {
        std::copy(to_copy.mp_alpha, to_copy.mp_alpha+m_ns, mp_alpha);
    }
    
    /**
     * Assignment operator.
     */
    ThirdbodyReactionStoich& operator = (ThirdbodyReactionStoich to_copy) {
        swap(*this, to_copy);
        return *this;
    }
    
    /**
     * Destructor.
     */
    ~ThirdbodyReactionStoich() {
        delete [] mp_alpha;
        mp_alpha = NULL;
    }
    
    /**
     * Adds this reaction's contribution to the species Jacobian matrix.  Note
     * that this will affect the rows of each species that particpates in this
     * reaction.  The compiler should be able to unroll each loop since the
     * sizes are known at compile time given the Reactants and Products 
     * stoichiometry types.
     */
    void contributeToJacobian(
        const double kf, const double kb, const double* const conc, 
        double* const work, double* const sjac, const size_t ns) const;
    
    friend void swap<Reactants, Products>(
        ThirdbodyReactionStoich<Reactants, Products>&,
        ThirdbodyReactionStoich<Reactants, Products>&);

private:

    size_t  m_ns;
    double* mp_alpha;
     
};


/**
 * This class is responsible for putting all of the pieces together and actually
 * computing the species jacobian.
 */
class JacobianManager
{
public:

    /**
     * Constructor.
     */
    JacobianManager(const Mutation::Thermodynamics::Thermodynamics& thermo)
        : m_thermo(thermo)
    {
        mp_work = new double [m_thermo.nSpecies()];
    }
    
    /**
     * Destructor.
     */
    ~JacobianManager() 
    {
        delete [] mp_work;
        
        std::vector<ReactionStoichBase*>::iterator iter = m_reactions.begin();
        for ( ; iter != m_reactions.end(); ++iter)
            delete *iter;
    }

    /**
     * Adds a reaction to the JacobianManager so that the contributions of the 
     * reaction will be included when computing the species Jacobian.
     */
    void addReaction(const Reaction& reaction);
    
    /**
     * Computes the square species source Jacobian matrix.
     */
    void computeJacobian(
        const double* const kf, const double* const kb, 
        const double* const conc, double* const sjac) const;
    
private:
    
    /**
     * Adds a new ReactionStoich object to the reaction list.  The reactant 
     * JacStoich** type is predetermined and input as a template parameter so 
     * that only the product type needs to be enumerated.
     */
    template <typename Reactants>
    void addReactionStoich(
        JacStoichBase* p_reacs, JacStoichBase* p_prods, const StoichType type, 
        const Reaction& reaction);

    /**
     * Determines which type of stoichiometry is defined by the set of species
     * names belonging to a reaction side and creates a new JaStoich** object
     * of that type.  Returns true if the stoichiometry is a valid type, false
     * otherwise.
     */
    bool getJacStoich(
        const std::vector<int>& stoich_set, JacStoichBase** p_stoich,
        StoichType& type) const;
    
private:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    double* mp_work;    
    std::vector<ReactionStoichBase*> m_reactions;

};

    } // namespace Kinetics
} // namespace Mutation

#endif // JACOBIAN_MANAGER_H



