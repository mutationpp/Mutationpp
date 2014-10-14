#include "Thermodynamics.h"
#include "MillikanWhite.h"
#include "TransferModel.h"
#include <cmath>

namespace Mutation {
    namespace Transfer {
      
// Implementation of the Vibrational-Translational Energy Transfer.
      
inline double const OmegaVT::compute_Park_correction_VT(int const i_vibrator, int const i_partner)
{
    // Limiting cross section for Park's Correction
    double m_transfer_P = mp_thermo->P(); 
    double m_transfer_T = mp_thermo->T();

    double sigma;
    if (m_transfer_T > 20000.0) {
      sigma = m_mw[i_vibrator].omega() * 6.25 ; // 6.25 = (50000/20000)^2
    } else {
      sigma = m_mw[i_vibrator].omega() *(2.5E9/(m_transfer_T*m_transfer_T));
    }
    
    return(m_const_Park_correction * sqrt(m_mw[i_vibrator][i_partner].mu()*m_transfer_T)/(sigma*m_transfer_P));
}
 
inline double const OmegaVT::compute_tau_VT_mj(int const i_vibrator, int const i_partner)
{
//    Enable in the future for multiple vibrational temperatures
      double m_transfer_P = mp_thermo->P();
      double m_transfer_T = mp_thermo->T();
  
      return( exp( m_mw[i_vibrator][i_partner].a() * (pow(m_transfer_T,  -1.E0/3.E0) - m_mw[i_vibrator][i_partner].b()) -18.42E0) * 101325.E0/ m_transfer_P );                                
}
      
double OmegaVT::compute_tau_VT_m(int const i_vibrator)
{
    const double * p_transfer_Y = mp_thermo->Y();
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    // Partner offset
    
        for (int i_partner = m_transfer_offset; i_partner < m_transfer_nHeavy; ++i_partner){
            double tau_j = compute_tau_VT_mj(i_vibrator, i_partner-m_transfer_offset) + compute_Park_correction_VT(i_vibrator, i_partner-m_transfer_offset);
            sum1 += p_transfer_Y[i_partner]/(mp_Mw[i_partner]);
            sum2 += p_transfer_Y[i_partner]/(mp_Mw[i_partner]*tau_j);
        }
  
    return(sum1/sum2);
    
}

          } // namespace Kinetics
} // namespace Mutation
      
      
      
      
      
      
      
      
      
      
