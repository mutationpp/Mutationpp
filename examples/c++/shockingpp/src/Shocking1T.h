#ifndef SHOCKING1T_H
#define SHOCKING1T_H

#include <boost/numeric/odeint.hpp>
#include "mutation++.h"

#include "Problem.h"

class Shocking1T: public Problem {
public:
    Shocking1T(Mutation::Mixture& l_mix, Data& l_data);
    ~Shocking1T(){}

    vector_type computedxdt(const vector_type& x, const double t);

private:
    size_t m_ns, m_nmom, m_nen, m_neq;
    size_t pos_u, pos_T;

    double m_mdot;

    mutable std::vector<double> v_Mw;
    mutable std::vector<double> v_yi, v_rhoi, v_omegai, v_hi, v_cpi;
    mutable vector_type v_dxdt;

    const size_t set_state_rhoi_T;
    const double m_energy_transfer_source;
};

#endif /* SHOCKING1T_H */
