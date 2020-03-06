/**
 * @file OmegaVT.cpp
 *
 * @brief Implementation of OmegaVT.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "MillikanWhite.h"
#include "Mixture.h"
#include "TransferModel.h"
#include <cmath>

using namespace Mutation;

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between vibrational and translational energy modes.
 */
class OmegaVT : public TransferModel
{
public:

    OmegaVT(Mixture& mix)
        : TransferModel(mix), m_mw(mix)
    {
        m_const_Park_correction = std::sqrt(PI*KB/(8.E0*NA));
        m_ns              = m_mixture.nSpecies();
        m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

        mp_Mw = new double [m_ns];
        for(int i = 0; i < m_ns; ++i)
            mp_Mw[i] = m_mixture.speciesMw(i);
        mp_hv = new double [m_ns];
        mp_hveq = new double [m_ns];
    }

    virtual ~OmegaVT()
    {
        delete [] mp_Mw;
        delete [] mp_hv;
        delete [] mp_hveq;
    }

    /**
     * Computes the source terms of the Vibration-Translational energy transfer in \f$ [J/(m^3\cdot s)] \f$
     * using a Landau-Teller formula taking into account Park's correction (default; can be disabled by making zero the appropriate flag, see below):
     * \f[ \Omega^{VT}_m = \rho_m \frac{e^V_m\left(T\right)-e^V_m\left(T_{vm}\right)} {\tau^{VT}_m} \f]
     *
     * The average relaxation time \f$ \tau^{VT}_m \f$ is given by the expression:
     *
     * \f[ \tau^{VT}_m = \frac{ \sum_{j \in \mathcal{H}} \rho_j / M_j}{ \sum_{j \in \mathcal{H}} \rho_j / (M_j \tau^{VT}_{mj}) } \f]
     *
     * More information about the above model can be found in @cite Panesi 2008
     *
     * @todo compare the average relaxation time with the one suggested by Park93
     *
     */

    double source()
    {
        const double * p_Y = m_mixture.Y();
        double rho = m_mixture.density();
        double T = m_mixture.T();
        double Tv = m_mixture.Tv();

        m_mixture.speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, mp_hveq, NULL, NULL);
        m_mixture.speciesHOverRT(T, Tv, T, Tv, Tv, NULL, NULL, NULL, mp_hv, NULL, NULL);

        int inv = 0;
        double src = 0.0;
        for (int iv = 0; iv-inv < m_mw.nVibrators(); ++iv){
            if(m_mixture.species(iv).type() != Mutation::Thermodynamics::MOLECULE){
                inv++;
            } else {
                src += p_Y[iv]*rho*RU*T/mp_Mw[iv]*(mp_hveq[iv] - mp_hv[iv])/compute_tau_VT_m(iv-inv);
            }
        }
        return src;
    }

private:

    MillikanWhite m_mw;

   /**
    * @brief This function computes the Millikan and White relaxation time
    * for each diatomic collision
    *
    * @param vibrator index
    * @param partner index
    *
    * @return Millikan and White relaxation time \tau_{m,j}^{MW}
    */

   inline double const compute_tau_VT_mj(int const, int const);

    /**
     * @brief Computes the frequency average over heavy particles.
     * @param vibrator index
     *
     * @return
     */
    double compute_tau_VT_m(int const);

    /**
     * @brief This function computes the Park correction
     * for vibrational-translational energy transfer
     * for each collision pair.
     *
     * @param vibrator index
     * @param partner index
     *
     * @return Park correction \tau_{m,j}^P
     */

    inline double const compute_Park_correction_VT(int const, int const);

    /**
     * Necessary variables
     */
    int m_ns;
    int m_transfer_offset;
    double* mp_Mw;
    double* mp_hv;
    double* mp_hveq;

    double m_const_Park_correction;
};
      
// Implementation of the Vibrational-Translational Energy Transfer.
      
inline double const OmegaVT::compute_Park_correction_VT(int const i_vibrator, int const i_partner)
{
    // Limiting cross section for Park's Correction
    double P = m_mixture.P();
    double T = m_mixture.T();

    double sigma;
    if (T > 20000.0) {
      sigma = m_mw[i_vibrator].omega() * 6.25 ; // 6.25 = (50000/20000)^2
    } else {
      sigma = m_mw[i_vibrator].omega() *(2.5E9/(T*T));
    }
    
    return(m_const_Park_correction * sqrt(m_mw[i_vibrator][i_partner].mu()*T)/(sigma*P));
}

 
inline double const OmegaVT::compute_tau_VT_mj(int const i_vibrator, int const i_partner)
{
//    Enable in the future for multiple vibrational temperatures
      double P = m_mixture.P();
      double T = m_mixture.T();

      return( exp( m_mw[i_vibrator][i_partner].a() * (pow(T,-1.0/3.0) - m_mw[i_vibrator][i_partner].b()) -18.421) * ONEATM / P );
}
      
double OmegaVT::compute_tau_VT_m(int const i_vibrator)
{
    const double * p_Y = m_mixture.Y();
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    // Partner offset
    for (int i_partner = m_transfer_offset; i_partner < m_ns; ++i_partner){
        double tau_j = compute_tau_VT_mj(i_vibrator, i_partner-m_transfer_offset) + compute_Park_correction_VT(i_vibrator, i_partner-m_transfer_offset);
        sum1 += p_Y[i_partner]/(mp_Mw[i_partner]);
        sum2 += p_Y[i_partner]/(mp_Mw[i_partner]*tau_j);
    }
  
    return(sum1/sum2);
}


// Register the transfer model
Utilities::Config::ObjectProvider<
    OmegaVT, TransferModel> omegaVT("OmegaVT");

    } // namespace Transfer
} // namespace Mutation 
