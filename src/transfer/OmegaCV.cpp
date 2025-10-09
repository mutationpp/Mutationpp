/**
 * @file OmegaCV.cpp
 *
 * @brief Implementation of OmegaCV.
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

#include "Mixture.h"
#include "TransferModel.h"

namespace Mutation {
    namespace Transfer {

/**
* Represents a coupling between chemistry and vibrational energy modes.
*/
class OmegaCV : public TransferModel
{
/**
* The necessary functions for the computation of the  source terms for
* Vibration-Chemistry-Vibration coupling are implemented.
*
* The dissociation of molecules is a procedure which removes energy from the vibrational mode
* of them and re-introduces it as chemical energy of the produced atoms and molecules. In order to accurately determine
* the exact energy lost from the vibrational mode in a dissociation reaction one must use state-to-state models in order to
* account for the fact that not all the vibrational states have the exact same probability to dissociate. Generally, the higher
* vibrational levels have a greater probability to dissociate compared to the lower ones because the latter ones should ladder climb
* before dissociating. Nonetheless, the state-to-state models are not easy to use in a CFD code due to the increase of the computational
* cost. Therefore, one should rely on accurate Boltzmann equilibrium models in order to take into account for this mode of energy transfer.
*
* Taking into account the above, the models for the Vibration-Chemistry-Vibration energy transfer can be generaly placed in to two different categories,
* the non-preferential and the preferential models. In the first category, all the energy states have the same probability to dissociate compared to the
* preferential ones which assume higher probability of dissociation for the higher energy levels. Currently, only the simplest non-preferential model has been
* implemented, but the code gives the opportunity to be extended for more complicated Vibration-Chemistry-Vibration energy transfer models.
*
*/

public:
    OmegaCV(Mutation::Mixture& mix)
        : TransferModel(mix)
    {
        m_ns = mixture().nSpecies();
        m_nt = mixture().nEnergyEqns();
        mp_wrk1 = new double [m_ns];
        mp_wrk2 = new double [m_ns];
        mp_wrk3 = new double [m_ns*m_ns];
        mp_wrk4 = new double [m_ns];
        mp_wrk5 = new double [m_ns];
        mp_wrk6 = new double [m_ns];
	};

    ~OmegaCV()
    {
        delete [] mp_wrk1;
        delete [] mp_wrk2;
        delete [] mp_wrk3;
        delete [] mp_wrk4;
        delete [] mp_wrk5;
        delete [] mp_wrk6;
    };
/**
 * Computes the source terms of the Vibration-Chemistry energy transfer in \f$ [J/(m^3\cdot s)] \f$
 *
 * \f[ \Omega^{CV}_{mi} = \sum c_1 e^V_{mi} \dot{\omega}_i\f]
 *
 */
    double source()
    {
        static int i_transfer_model = 0;
        switch (i_transfer_model){
            case 0:
              return compute_source_Candler();
            break;
            default:
                std::cerr << "The selected Chemistry-Vibration-Chemistry model is not implemented yet";
                return 0.0;
        }
    }

/**
 * Computes the density-based jacobian terms of the Vibration-Chemistry energy transfer.$
 *
 * \f[ \frac{\partial \Omega^{CV}_{mi}}{\partial \rho_{i}} = \sum c_1 e^V_{mi} \frac{\partial \dot{\omega}_i}{\partial \rho_{i}} \f]
 *
 */
    void jacobianRho(double* const p_jacRho);

/**
 * Computes the temperature-based jacobian terms of the Vibration-Chemistry energy transfer.$
 *
 * \f[ \frac{\partial \Omega^{CV}_{mi}}{\partial T} = \sum c_1 \frac{\partial e^V_{mi}}{\partial T} \dot{\omega}_i 
 *                                                     + \sum c_1 e^V_{mi} \frac{\partial \dot{\omega}_i}{\partial T} \f]
 */
    void jacobianTTv(double* const p_jacTTv);

private:
    int m_ns;
    int m_nt;
    double* mp_wrk1;
    double* mp_wrk2;
    double* mp_wrk3;
    double* mp_wrk4;
    double* mp_wrk5;
    double* mp_wrk6;

    double const compute_source_Candler();
};

 /**
 * Non-preferential Model according to Candler with
 *
 * with \f$ c_1 \f$ equal to 1 for non-preferential models. For more information:
 *
 * G. V. Candler, MacCormack, Computation of weakly ionized hypersonic flows in thermochemical nonequilibrium, Journal of
 * Thermophysics and Heat Transfer, 1991, 5(11):266
 *
 */

    double const OmegaCV::compute_source_Candler()
    {
        
        
        // Getting Vibrational Energy
        mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);
        
        // Getting Production Rate
        mixture().netProductionRates(mp_wrk2);
        
        // Inner Product
        double c1 = 1.0E0;
        double sum = 0.E0;
        
        for(int i = 0 ; i < m_ns; ++i)
            sum += mp_wrk1[i]*mp_wrk2[i]/mixture().speciesMw(i);
        
        return(c1*sum*mixture().T()*RU);
    }

/**
 * Assuming non-preferential model for now, density-based jacobians are computed.
 *
 */

    void OmegaCV::jacobianRho(double* const p_jacRho)
    {
        std::fill(p_jacRho, p_jacRho + m_ns, 0.);
        
        // Getting Vibrational Energy
        mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);
        
        // Getting density based jacobian for production rates
        mixture().jacobianRho(mp_wrk3);
        
        //Inner product
        double c1 = 1.0E0;
        double aux = c1*RU*mixture().T();
        
        for(int i = 0; i < m_ns; ++i)
            for(int j = 0; j < m_ns; ++j)
                p_jacRho[i] += aux*mp_wrk1[i]*mp_wrk3[j*m_ns+i]/mixture().speciesMw(i);
    }

/**
 * Assuming non-preferential model for now, temperature-based jacobians are computed.
 *
 */

    void OmegaCV::jacobianTTv(double* const p_jacTTv)
    {
        std::fill(p_jacTTv, p_jacTTv + m_nt, 0.);
        
        double T=mixture().T();
        double Tv=mixture().Tv();
        
        // Getting Vibrational Energy
        mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);
        
        // Getting production rates
        mixture().netProductionRates(mp_wrk2);
        
        // Getting vibrational energy derivative with respect to Tv
        mixture().speciesCpOverR(T, Tv, T, Tv, Tv, NULL, NULL, NULL, mp_wrk4, NULL);
        
        // Getting temperature based jacobians for production rates
        mixture().jacobianT(mp_wrk5);
        mixture().jacobianTv(mp_wrk6);
        
        //Inner product
        double c1 = 1.0E0;
        double aux = c1*RU*T;
        double aux1 = c1*RU;
        
        for(int i = 0; i < m_ns; ++i) {
            p_jacTTv[0] += aux*mp_wrk1[i]*mp_wrk5[i]/mixture().speciesMw(i);
            p_jacTTv[1] += aux*mp_wrk1[i]*mp_wrk6[i]/mixture().speciesMw(i);
            p_jacTTv[1] += aux1*mp_wrk4[i]*mp_wrk2[i]/mixture().speciesMw(i);
        }
    }

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCV, TransferModel> omegaCV("OmegaCV");
      
    } // namespace Transfer
} // namespace Mutation
