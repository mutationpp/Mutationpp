#include "LarsenTTv.h"

LarsenTTv::LarsenTTv(Mutation::Mixture& l_mix, Data& l_data)
            : Problem(l_mix, l_data),
              m_ns(l_mix.nSpecies()),
              m_nen(l_mix.nEnergyEqns()),
              m_neq(m_ns+m_nen),
              v_Mw(m_ns, 0.),
              v_yi(m_ns, 0.),
              v_rhoi(m_ns, 0.),
              v_omegai(m_ns, 0.),
              v_T(m_nen, 0.),
              v_hi(m_ns*m_nen, 0.),
              v_ei(m_ns*m_nen, 0.),
              v_cpi(m_ns*m_nen, 0.),
              v_cvi(m_ns*m_nen, 0.),
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1),
              v_omega_int_en(m_nen-1,0.),
              m_energy_transfer_source(0.0){
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

vector_type LarsenTTv::computedxdt(const vector_type& xState, const double t){

  double x_ext[2];
  double u_ext[2];
  double rho_ext[2];
  double dummy[2];

  double u_now, rho_now;

  // Initializing variables to current state
  m_data.getStepValues(x_ext, u_ext, rho_ext);

  u_now   = u_ext[0] + (u_ext[1] - u_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);
  rho_now = rho_ext[0] + (rho_ext[1] - rho_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);


  // Compute densities..
  for(size_t i_sp = 0; i_sp < m_ns; ++i_sp) {
    v_yi[i_sp] = xState[i_sp];
    if(xState[i_sp] < 0.0){   // It might become negative during some kind of integration..
      v_yi[i_sp] = 0.0;
    }
    v_rhoi[i_sp] = v_yi[i_sp] * rho_now;
  }

  // right after the m_ns mass equations, there are temperatures
  for(size_t i_en = 0; i_en < m_nen; ++i_en) {
    v_T[i_en] = xState[m_ns + i_en];
  }

  // ----  Set mixture state and get properties from the mixture  ----
  // NOTE: &rhoi[0] is a tricky way to pass a VECTOR where an ARRAY was asked..
  m_mix.setState(&v_rhoi[0], &v_T[0], set_state_rhoi_T);

  m_mix.getEnthalpiesMass(&v_hi[0]);
  m_mix.getEnergiesMass(&v_ei[0]);
  m_mix.getCpsMass(&v_cpi[0]);
  m_mix.getCvsMass(&v_cvi[0]);
  m_mix.netProductionRates(&v_omegai[0]);
  m_mix.energyTransferSource(&v_omega_int_en[0]);

  // ----  Computing derivatives  -----
  // ====> Mass equations
  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_dxdt[ii] = v_omegai[ii]/rho_now;
  }

  // ====> Vibrational Temperature
  double fac_num = 0.0;
  double fac_den = 0.0;
  for (int i_sp = 0; i_sp < m_ns; ++i_sp ){
      fac_num += v_omegai[i_sp] * v_ei[i_sp + m_ns];
      fac_den += v_yi[i_sp] * v_cvi[i_sp + m_ns];
  }

  v_dxdt[m_ns+1] = (v_omega_int_en[0] - fac_num) / (fac_den * rho_now);

  // ====> Translational Temperature
  double num_part1 = 0.0;
  double num_part2 = 0.0;
  double den       = 0.0;
  for(size_t ii = 0; ii < m_ns; ++ii) {
    num_part1 += v_hi[ii]*v_omegai[ii]/rho_now;
    num_part2 += v_yi[ii]*v_cpi[ii+m_ns]*v_dxdt[m_ns+1];
    den       += v_yi[ii]*v_cpi[ii];
  }
  // Approximate derivative of u^2
  double Du2Dt = u_now*(pow(u_ext[1],2) - pow(u_ext[0],2))/(x_ext[1] - x_ext[0] + 1E-15);

  v_dxdt[m_ns] = - (0.5*Du2Dt + num_part1 + num_part2)/den;

  // Correcting ALL THE derivatives with the velocity (remember, these are material derivatives!!)
  for(size_t ii = 0; ii < m_ns+m_nen; ++ii) {
    v_dxdt[ii] /= u_now;
  }

  return v_dxdt;

};
                    
