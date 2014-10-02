#include "Thermodynamics.h"
#include "MillikanWhite.h"
#include "TransferModel.h"
#include <cmath>

namespace Mutation {
    namespace Transfer {
      
// Implementation of the Vibrational-Translational Energy Transfer.
      
inline double const TransferModelVT::compute_Park_correction_VT(int const i_vibrator, int const i_partner)
{
    // Limiting cross section for Park's Correction
    double sigma = m_mw[i_vibrator].omega() *(2.5E9/(m_transfer_T*m_transfer_T)); 
    return(m_const_Park_correction * sqrt(m_mw[i_vibrator][i_partner].mu())/sigma);
}
 
inline double const TransferModelVT::compute_tau_VT_mj(int const i_vibrator, int const i_partner)
{
//    Enable in the future for multiple vibrational temperatures
      return( exp( m_mw[i_vibrator][i_partner].a() * (pow(m_transfer_T,  -1.E0/3.E0)) - m_mw[i_vibrator][i_partner].b() -18.42E0) * 101325.E0/ m_transfer_P );                                
}
      
double TransferModelVT::compute_tau_VT_m()
{
        
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    
    // Partner offset
    
    for (int i_vibrator = 0; i_vibrator < m_mw.nVibrators(); ++i_vibrator){
        for (int i_partner = 0; i_partner < m_transfer_nHeavy; ++i_partner){
            double tau_j = compute_tau_VT_mj(i_vibrator, i_partner) + compute_Park_correction_VT(i_vibrator, i_partner);
             
            sum1 += p_transfer_Y[i_partner]/m_mw[i_vibrator][i_partner].mu();
            sum2 += p_transfer_Y[i_partner]/(m_mw[i_vibrator][i_partner].mu() * tau_j);
            sum3 += p_transfer_Y[i_partner];
            
        }
    }
  
    double tau_m = sum1/sum2;
  
    return(tau_m);
    
}
      
          } // namespace Kinetics
} // namespace Mutation
      
      
      
      
      
      
      
      
      
      