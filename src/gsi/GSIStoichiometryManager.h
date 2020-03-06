/**
 * @file GSIStoichiometryManager.h
 *
 * @brief Declaration of GSIStoichiometryManager class.
 */

/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
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


#ifndef GSI_STOICHIOMETRY_MANAGER_H
#define GSI_STOICHIOMETRY_MANAGER_H

#include <Eigen/Dense>

namespace Mutation {
    namespace GasSurfaceInteraction {

template <int N>
class Stoich
{
public:
    Stoich(){
        m_rxn = 0;
        for (int i = 0; i < N; ++i) m_sps[i] = 0;
    }

    inline void multReaction(
        const Eigen::VectorXd& p_s, Eigen::VectorXd& p_r) const
    {
        for (int i = 0; i < N; ++i) p_r(m_rxn) *= p_s(m_sps[i]);
    }

    inline void incrReaction(
        const Eigen::VectorXd& p_s, Eigen::VectorXd& p_r) const
    {
        for (int i = 0; i < N; ++i) p_r(m_rxn) += p_s(m_sps[i]);
    }

    inline void decrReaction(
        const Eigen::VectorXd& p_s, Eigen::VectorXd& p_r) const
    {
        for (int i = 0; i < N; ++i) p_r(m_rxn) -= p_s(m_sps[i]);
    }

    inline void incrSpecies(
        const Eigen::VectorXd& p_r, Eigen::VectorXd& p_s) const
    {
        for (int i = 0; i < N; ++i) p_s(m_sps[i]) += p_r(m_rxn);
    }

    inline void decrSpecies(
        const Eigen::VectorXd& p_r, Eigen::VectorXd& p_s) const
    {
        for (int i = 0; i < N; ++i) p_s(m_sps[i]) -= p_r(m_rxn);
    }

protected:
    size_t m_rxn;
    size_t m_sps[N];

}; //class Stoich<N>

//==============================================================================

class Stoich1 : public Stoich<1>
{
public:
    Stoich1(size_t rxn, size_t sp1) {
        m_rxn = rxn; m_sps[0] = sp1;
    }
};

//==============================================================================

class Stoich2 : public Stoich<2>
{
public:
    Stoich2(size_t rxn, size_t sp1, size_t sp2)
    {
        m_rxn = rxn;
        if (sp1 <= sp2) {
            m_sps[0] = sp1; m_sps[1] = sp2;
        } else {
            m_sps[0] = sp2; m_sps[1] = sp1;
        }
    }
};

//==============================================================================

class Stoich3 : public Stoich<3>
{
public:
    Stoich3(size_t rxn, size_t sp1, size_t sp2, size_t sp3)
    {
        m_rxn = rxn;
        if (sp1 <= sp2) {
            m_sps[0] = sp1; m_sps[1] = sp2;
        } else {
            m_sps[0] = sp2; m_sps[1] = sp1;
        }

        if (m_sps[1] <= sp3) m_sps[2] = sp3;
        else {
            if (sp3 <= m_sps[0]) {
                m_sps[2] = m_sps[1]; m_sps[1] = m_sps[0]; m_sps[0] = sp3;
            } else {
                m_sps[2] = m_sps[1]; m_sps[1] = sp3;
            }
        }
    }
};


//==============================================================================

/**
 * Provides sparse-matrix type operations for quantities that depend
 * on reaction stoichiometries.
 */
class GSIStoichiometryManager{
public:
    GSIStoichiometryManager( ){ }

    void addReaction(const int rxn, const std::vector<int>& sps);
    virtual void multReactions(
        const Eigen::VectorXd& p_s, Eigen::VectorXd& p_r) const;
    virtual void incrReactions(
        const Eigen::VectorXd& p_s, Eigen::VectorXd& p_r) const;
    virtual void decrReactions(
        const Eigen::VectorXd& p_s, Eigen::VectorXd& p_r) const;
    virtual void incrSpecies(
        const Eigen::VectorXd& p_r, Eigen::VectorXd& p_s) const;
    virtual void decrSpecies(
        const Eigen::VectorXd& p_r, Eigen::VectorXd& p_s) const;

    virtual ~GSIStoichiometryManager(){ }

private:
    std::vector<Stoich1> m_stoich1_vec;
    std::vector<Stoich2> m_stoich2_vec;
    std::vector<Stoich3> m_stoich3_vec;

};

    }// namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSI_STOICHIOMETRY_MANAGER_H
