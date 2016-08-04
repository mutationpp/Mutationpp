#ifndef LARSEN1T_H
#define LARSEN1T_H

#include <boost/numeric/odeint.hpp>
#include "mutation++.h"

#include "Problem.h"

class Larsen1T: public Problem {
public:
    Larsen1T(Mutation::Mixture& l_mix, Data& l_data);
    ~Larsen1T(){}

    vector_type computedxdt(const vector_type& x, const double t);

private:
    size_t m_ns, m_nen, m_neq;

    mutable std::vector<double> v_Mw;
    mutable std::vector<double> v_yi, v_rhoi, v_omegai, v_hi, v_cpi;
    mutable vector_type v_dxdt;

    const size_t set_state_rhoi_T;
};

#endif /* LARSEN1T_H */
