#ifndef LARSENTTv_H
#define LARSENTTv_H

#include <boost/numeric/odeint.hpp>
#include "mutation++.h"

#include "Problem.h"

class LarsenTTv: public Problem {
public:
    LarsenTTv(Mutation::Mixture& l_mix, Data& l_data);
    ~LarsenTTv(){}

    vector_type computedxdt(const vector_type& x, const double t);

private:
    size_t m_ns, m_nen, m_neq;

    mutable std::vector<double> v_Mw;
    mutable std::vector<double> v_yi, v_rhoi, v_omegai, v_hi, v_cpi, v_ei, v_cvi, v_omega_int_en, v_T;
    mutable vector_type v_dxdt;

    const size_t set_state_rhoi_T;
    const double m_energy_transfer_source;
};

#endif /* LARSENTTv_H */
