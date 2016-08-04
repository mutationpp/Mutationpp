#include "Larsen1T.h"

Larsen1T::Larsen1T(Mutation::Mixture& l_mix, Data& l_data)
            : Problem(l_mix, l_data),
              m_ns(l_mix.nSpecies()),
              m_nen(l_mix.nEnergyEqns()),
              m_neq(m_ns+m_nen),
              v_Mw(m_ns, 0.),
              v_yi(m_ns, 0.),
              v_rhoi(m_ns, 0.),
              v_omegai(m_ns, 0.),
              v_hi(m_ns, 0.),
              v_cpi(m_ns, 0.),
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1){
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

vector_type Larsen1T::computedxdt(const vector_type& xState, const double t){

  double x_ext[2];
  double u_ext[2];
  double rho_ext[2];

  double u_now, rho_now, T_now;

  // Initializing variables to current state
  m_data.getStepValues(x_ext, u_ext, rho_ext);

  u_now   = u_ext[0] + (u_ext[1] - u_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);
  rho_now = rho_ext[0] + (rho_ext[1] - rho_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);


  // Compute densities..
  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_yi[ii] = xState[ii];
    if(xState[ii] < 0.0){   // It might become negative during some kind of integration..
      v_yi[ii] = 0.0;
    }
    v_rhoi[ii] = v_yi[ii] * rho_now;
  }

  T_now = xState[m_ns];

  // Set mixture state
  // NOTE: &rhoi[0] is a tricky way to pass a VECTOR where an ARRAY was asked..
  m_mix.setState(&v_rhoi[0], &T_now, 1);

  // ----  Get properties from the mixture  ----
  m_mix.netProductionRates(&v_omegai[0]);
  m_mix.speciesHOverRT(T_now, &v_hi[0]);
  m_mix.speciesCpOverR(T_now, &v_cpi[0]);

  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_hi[ii]  *= Mutation::RU*T_now/m_mix.speciesMw(ii);
    v_cpi[ii] *= Mutation::RU/m_mix.speciesMw(ii);
  }

  // ----  Computing derivatives  -----
  // ====> Mass equations
  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_dxdt[ii] = v_omegai[ii]/rho_now;
  }

  // ====> Translational temperature
  double num_part = 0.0;
  double den      = 0.0;
  for(size_t ii = 0; ii < m_ns; ++ii) {
    num_part += v_hi[ii]*v_omegai[ii]/rho_now;
    den      += v_yi[ii]*v_cpi[ii];
  }
  // Approximate derivative of u^2
  double Du2Dt = u_now*(pow(u_ext[1],2) - pow(u_ext[0],2))/(x_ext[1] - x_ext[0] + 1E-15);

  v_dxdt[m_ns] = - (0.5*Du2Dt + num_part)/den;

  // Correcting ALL THE derivatives with the velocity (remember, these are material derivatives!!)
  for(size_t ii = 0; ii < m_ns+ 1; ++ii) {
    v_dxdt[ii] /= u_now;
  }

  return v_dxdt;

};
                    
