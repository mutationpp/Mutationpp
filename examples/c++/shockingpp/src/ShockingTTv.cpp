#include "ShockingTTv.h"

ShockingTTv::ShockingTTv(Mutation::Mixture& l_mix, Data& l_data)
            : Problem(l_mix, l_data),
              m_ns(l_mix.nSpecies()),
              m_nmom(1),
              m_nen(l_mix.nEnergyEqns()),
              m_neq(m_ns+m_nmom+m_nen),
              pos_u(m_ns),
              pos_T(m_ns+m_nmom),
              pos_T_int(pos_T+1),
              m_mdot(0.),
              v_Mw(m_ns, 0.),
              v_yi(m_ns, 0.),
              v_rhoi(m_ns, 0.),
              v_omegai(m_ns, 0.),
              v_omega_int_en(m_nen-1, 0.),
              v_T(m_nen, 0.),
              v_hi(m_ns*m_nen, 0.),
              v_ei(m_ns*m_nen, 0.),
              v_cpi(m_ns*m_nen, 0.),
              v_cvi(m_ns*m_nen, 0.),
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1),
              m_energy_transfer_source(0.0){
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

vector_type ShockingTTv::computedxdt(const vector_type& x, const double t){

       m_mdot = m_data.mdot();
       double l_u = x[pos_u];
       double l_rho = m_mdot/l_u;

       for (int i_sp = 0; i_sp < m_ns; i_sp++){
           v_yi[i_sp] = x[i_sp];
           v_rhoi[i_sp] = v_yi[i_sp] * l_rho;
       }
       for (int i_en = 0; i_en < m_nen; ++i_en){
           v_T[i_en] = x[pos_T+i_en];
       }

    m_mix.setState(&v_rhoi[0], &v_T[0], set_state_rhoi_T);
    m_mix.getEnthalpiesMass(&v_hi[0]);
    m_mix.getEnergiesMass(&v_ei[0]);
    m_mix.getCpsMass(&v_cpi[0]);
    m_mix.getCvsMass(&v_cvi[0]);
    m_mix.netProductionRates(&v_omegai[0]);
    m_mix.energyTransferSource(&v_omega_int_en[0]);

    // Sum Variables
    double wrk1, wrk2, wrk3;
    double h_sum =0.0;
    double cp_sum = 0.0;
    double y_sum = 0.0;
    double omega_sum = 0.0;

    for (int i_sp = 0; i_sp < m_ns; i_sp++){
        wrk1 = 1.0/v_Mw[i_sp];
        wrk2 = v_yi[i_sp];
        wrk3 = v_omegai[i_sp];

        h_sum += v_hi[i_sp]*wrk3;
        cp_sum += wrk2 * v_cpi[i_sp];
        y_sum += wrk2 * wrk1;
        omega_sum += wrk3 * wrk1;
    }

    // Helper Variables;
    double a, b, c, d, e, f;
    a = m_mdot * l_u / (Mutation::RU * v_T[0]) - l_rho * y_sum;
    b = m_mdot / v_T[0] * y_sum;
    c = - omega_sum;
    d = m_mdot * l_u;
    e = m_mdot * cp_sum;
    f = -h_sum + m_mdot * m_energy_transfer_source; // This is for radiation

    for (int i_en = 0; i_en < m_nen-1; i_en++){
        f = f - v_omega_int_en[i_en];
    }

    wrk1 = 1.0/m_mdot;
    wrk2 = 1.0/(a*e - b*d);

    // Filling the derivatives vector
    for (int i_sp = 0; i_sp < m_ns; i_sp++){
        v_dxdt[i_sp] = v_omegai[i_sp] * wrk1;
    }

    v_dxdt[pos_u] = (c*e - b*f) * wrk2;
    v_dxdt[pos_T] = (a*f - c*d) * wrk2;

    // This is the part that need to be generalized in order to be more than TTv
    // Probably a function in Thermodynamics should be added that computes w*h
    // Fix it for electrons
    double fac_num = 0.;
    double fac_den = 0.;
    for (int i_sp = 0; i_sp < m_ns; ++i_sp ){
        fac_num += v_omegai[i_sp] * v_ei[i_sp + m_ns];
        fac_den += v_yi[i_sp] * v_cvi[i_sp + m_ns];
    }

    for (int i_en = 0; i_en < m_nen - 1; ++i_en){
        v_dxdt[pos_T_int + i_en] = (v_omega_int_en[i_en] - fac_num) / fac_den * wrk1;
    }

    return v_dxdt;
}
