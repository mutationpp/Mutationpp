#ifndef SHOCKINGNT_H
#define SHOCKINGNT_H

#include <boost/numeric/odeint.hpp>
#include "mutation++.h"

#include "Problem.h"

class ShockingTTv: public Problem {
public:
    ShockingTTv(Mutation::Mixture& l_mix, Data& l_data);
    ~ShockingTTv(){}

    vector_type computedxdt(const vector_type& x, const double t);

private:
    size_t m_ns, m_nmom, m_nen, m_neq;
    size_t pos_u, pos_T, pos_T_int;

    double m_mdot;

    mutable std::vector<double> v_Mw;
    mutable std::vector<double> v_yi, v_rhoi, v_omegai, v_T,  v_hi, v_cpi, v_omega_int_en;
    mutable std::vector<double> v_hti, v_hri, v_hvi, v_heli;
    mutable std::vector<double> v_cpti, v_cvri, v_cvvi, v_cveli;
    mutable vector_type v_dxdt;

    const size_t set_state_rhoi_T;
    const double m_energy_transfer_source;
};

#endif /* SHOCKINGNT_H */
