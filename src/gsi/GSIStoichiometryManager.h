#ifndef GSI_STOICHIOMETRY_MANAGER_H
#define GSI_STOICHIOMETRY_MANAGER_H

#include <vector>

#include "GSIStoichiometryManager.h"

namespace Mutation {
    namespace gsi {

//=======================================================================================

template <int N>
class Stoich
{
public:

    Stoich(){
        m_rxn = 0;
        for (int i = 0; i < N; ++i) m_sps[i] = 0;
    }

    inline void multReaction(const double* const p_s, double* const p_r) const {
        for (int i = 0; i < N; ++i) p_r[m_rxn] *= p_s[m_sps[i]];
    }
    inline void incrReaction(const double* const p_s, double* const p_r) const {
        for (int i = 0; i < N; ++i) p_r[m_rxn] += p_s[m_sps[i]];
    }
    inline void decrReaction(const double* const p_s, double* const p_r) const {
        for (int i = 0; i < N; ++i) p_r[m_rxn] -= p_s[m_sps[i]];
    }
    inline void incrSpecies(const double* const p_r, double* const p_s) const {
        for (int i = 0; i < N; ++i) p_s[m_sps[i]] += p_r[m_rxn];
    }
    inline void decrSpecies(const double* const p_r, double* const p_s) const {
        for (int i = 0; i < N; ++i) p_s[m_sps[i]] -= p_r[m_rxn];
    }

protected:

    size_t m_rxn;
    size_t m_sps[N];

}; //class Stoich<N>

//=======================================================================================

class Stoich1 : public Stoich<1> {
public:
    Stoich1(size_t rxn, size_t sp1) {
        m_rxn = rxn; m_sps[0] = sp1;
    }
};

//=======================================================================================

class Stoich2 : public Stoich<2> {
public:
    Stoich2(size_t rxn, size_t sp1, size_t sp2) {
        m_rxn = rxn;
        if (sp1 <= sp2) {
            m_sps[0] = sp1; m_sps[1] = sp2;
        } else {
            m_sps[0] = sp2; m_sps[1] = sp1;
        }
    }
};

//=======================================================================================

class Stoich3 : public Stoich<3>
{
public:
    Stoich3(size_t rxn, size_t sp1, size_t sp2, size_t sp3) {
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


//=======================================================================================

/**
 * Provides sparse-matrix type operations for quantities that depend on reaction
 * stoichiometries.
 */
class GSIStoichiometryManager{
public:
    GSIStoichiometryManager( ){ }

    virtual void addReaction(const int, const std::vector<int>&) = 0;
    virtual void multReactions(const double* const, double* const) const = 0;
    virtual void incrReactions(const double* const, double* const) const = 0;
    virtual void decrReactions(const double* const, double* const) const = 0;
    virtual void incrSpecies(const double* const, double* const) const = 0;
    virtual void decrSpecies(const double* const, double* const) const = 0;

    virtual ~GSIStoichiometryManager(){ }

};

//=======================================================================================

class GSIStoichiometryManager_GammaModel: public GSIStoichiometryManager{
public:
    GSIStoichiometryManager_GammaModel( ){ }

    void addReaction(const int rxn, const std::vector<int>& sps);
    void multReactions(const double* const p_s, double* const p_r) const;
    void incrReactions(const double* const p_s, double* const p_r) const;
    void decrReactions(const double* const p_s, double* const p_r) const;
    void incrSpecies(const double* const p_r, double* const p_s) const;
    void decrSpecies(const double* const p_r, double* const p_s) const;

    ~GSIStoichiometryManager_GammaModel(){ }

private:
    std::vector<Stoich1> m_stoich1_vec;
    std::vector<Stoich2> m_stoich2_vec;

};

//=======================================================================================

class GSIStoichiometryManager_FRC: public GSIStoichiometryManager{
public:
    GSIStoichiometryManager_FRC(){ }

    void addReaction(const int rxn, const std::vector<int>& sps);
    void multReactions(const double* const p_s, double* const p_r) const;
    void incrReactions(const double* const p_s, double* const p_r) const;
    void decrReactions(const double* const p_s, double* const p_r) const;
    void incrSpecies(const double* const p_r, double* const p_s) const;
    void decrSpecies(const double* const p_r, double* const p_s) const;

    ~GSIStoichiometryManager_FRC(){ }

private:
    std::vector<Stoich1> m_stoich1_vec;
    std::vector<Stoich2> m_stoich2_vec;
    std::vector<Stoich3> m_stoich3_vec;

    //std::vector<StoichJac#> m_stoich_jac1_vec;

};

//=======================================================================================

    }// namespace gsi
} // namespace Mutation

#endif // GSI_STOICHIOMETRY_MANAGER_H
