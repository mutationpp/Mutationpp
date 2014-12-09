/**
 * @file TransferModel.h
 *
 * @brief Energy Transfer Models
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef TRANSFER_TRANSFER_MODEL_H
#define TRANSFER_TRANSFER_MODEL_H

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "Transport.h"
#include "MillikanWhite.h"

namespace Mutation {
    namespace Transfer {

/**
 * @defgroup transfermodels Energy Transfer Models
 * @{
 */

class TransferModel
{
public:
    virtual ~TransferModel() { }
    //virtual void source(double* const p_source) = 0; 
    virtual double source() = 0; 
};

/**
 * Represents a coupling between vibrational and translational energy modes.
 */
class OmegaVT : public TransferModel
{
public:

    OmegaVT(const Thermodynamics::Thermodynamics& thermo)
        : m_mw(thermo){ 
            mp_thermo = &thermo;
            m_const_Park_correction = sqrt(PI* KB / (8.E0 * NA));

            m_ns              = thermo.nSpecies();
            m_transfer_nHeavy = thermo.nHeavy();
            m_transfer_offset = thermo.hasElectrons() ? 1 : 0;
 
	    mp_Mw = new double [m_ns];
            for(int i = 0; i < m_ns; ++i)
                mp_Mw[i] = thermo.speciesMw(i);               
            mp_hv = new double [m_ns];
            mp_hveq = new double [m_ns];
    }
 
    /**
     * Computes the source terms of the Vibration-Translational energy transfer in \f$ [J/(m^3\cdot s)] \f$
     * using a Landau-Teller formula taking into account Park's correction (default; can be disabled by making zero the appropriate flag, see below):
     * \f[ \Omega^{VT}_m = \rho_m \frac{e^V_m\left(T\right)-e^V_m\left(T_{vm}\right)}
     *      {\tau^{VT}_m} \left|\frac{T_{sh}-T_{V}}{T_{sh}-T_{Vsh}}\right|^{s-1} \f]
     * with \f$ s = 3.5exp(-\frac{5000}{T_s}) \f$. 
     * 
     * The average relaxation time \f$ \tau^{VT}_m \f$ is given by the expression:
     *
     * \f[ \tau^{VT}_m = \frac{ \sum_{j \in \mathcal{H}} \rho_j / M_j}{ \sum_{j \in \mathcal{H}} \rho_j / (M_j \tau^{VT}_{mj}) } \f]
     *
     * More information about the above model can be found in @cite Park1993 and 
     * @cite Schwarz1952 .
     * 
     */

    double source()
    { 
        const double * p_Y = mp_thermo->Y();
        double rho = mp_thermo->density();
        double T = mp_thermo->T();
        double Tv = mp_thermo->Tv();
    
        mp_thermo->speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, mp_hveq, NULL, NULL);
        mp_thermo->speciesHOverRT(T, Tv, T, Tv, Tv, NULL, NULL, NULL, mp_hv, NULL, NULL);

        int inv = 0;
        double src = 0.0;
        for (int iv = 0; iv-inv < m_mw.nVibrators(); ++iv){
            if(mp_thermo->species(iv).type() != Mutation::Thermodynamics::MOLECULE){
                inv++; 
	    } else {
                src += p_Y[iv]*rho*RU*T/mp_Mw[iv]*(mp_hveq[iv] - mp_hv[iv])/compute_tau_VT_m(iv-inv);
	    }
        }
        return src;
    }

private:
    const Thermodynamics::Thermodynamics* mp_thermo;
    MillikanWhite m_mw;
    
    // Declaring a function to compute \tau_VT_mj for each species. It returns a pointer to tau_VT_m.
    inline double const compute_tau_VT_mj(int const, int const);
    
    /**
     * Computes the frequency avarage over heavy particles.
     */
    double compute_tau_VT_m(int const);
    
    /**
     * This function computes the Park correction for vibrational-translational energy trasfer.
     */
    inline double const compute_Park_correction_VT(int const, int const);
    
    /**
     * Necesseary variables
     */
    int m_ns;
    int m_transfer_nHeavy; 
    int m_transfer_offset;
    double* mp_Mw;
    double* mp_hv;
    double* mp_hveq;
    
    double m_const_Park_correction;
};

/**
 * Represents a coupling between chemistry and vibational energy modes.
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
    OmegaCV(const Thermodynamics::Thermodynamics& thermo, Kinetics::Kinetics& kinetics){
        mp_thermo = &thermo; 
        mp_kinetics = &kinetics;
        m_ns = thermo.nSpecies();
        mp_wrk1 = new double [m_ns];
        mp_wrk2 = new double [m_ns];
    };

    ~OmegaCV()
    {
        delete [] mp_wrk1;
        delete [] mp_wrk2;
    };
      
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
    };
    
private:
    const Thermodynamics::Thermodynamics* mp_thermo;
    Kinetics::Kinetics* mp_kinetics;
    int m_ns;
    double* mp_wrk1;
    double* mp_wrk2;

    double const compute_source_Candler();
};

/**
 * Represents a coupling between chemistry and electronic energy modes.
 */
class OmegaCElec : public TransferModel
{
   /**
    * The necessary functions for the computation of the  source terms for 
    * Electronic-Chemistry-Electronic coupling are implemented. 
    * 
    */

public:
    OmegaCElec(const Thermodynamics::Thermodynamics& thermo, Kinetics::Kinetics& kinetics){
        mp_thermo = &thermo; 
        mp_kinetics = &kinetics;
        m_ns = thermo.nSpecies();
        mp_wrk1 = new double [m_ns];
        mp_wrk2 = new double [m_ns];
    };

    ~OmegaCElec()
    {
        delete [] mp_wrk1;
        delete [] mp_wrk2;
    };
      
    double source()
    {
      static int i_transfer_model = 0;
      switch (i_transfer_model){
          case 0:
              return compute_source_Candler();
	      break;
          default:
              std::cerr << "The selected Electronic-Chemistry-Electronic model is not implemented yet";
              return 0.0;
      }
    };
    
private:
    const Thermodynamics::Thermodynamics* mp_thermo;
    Kinetics::Kinetics* mp_kinetics;
    int m_ns;
    double* mp_wrk1;
    double* mp_wrk2;

    double const compute_source_Candler();
};

class OmegaET : public TransferModel
{
public:
    OmegaET(const Thermodynamics::Thermodynamics& thermo, Transport::Transport& transport){
        mp_thermo = &thermo;
        mp_transport = &transport;
        mp_collisions = transport.collisionData();
    };

    /**
     * Computes the source term of the Electron-Heavy Translational energy transfer in \f$ [J/(m^3\cdot s)] \f$
     * using a Landau-Teller formula:
     * \f[ \Omega^{ET} = \rho_e \frac{e^T_e\left(T\right)-e^T_e\left(T_{e}\right)}
     *      {\tau^{ET}}  \f]
     * 
     * The relaxation time \f$ \tau^{ET} \f$ is given by the expression:
     *
     * \f[ \tau^{ET} = \left( \sum_{j \in \mathcal{H}} \frac{2 M_e}{ M_i} \nu_{ei} \right)^{-1} \f]
     *
     */
    
    double source()
    {
        const double * p_Y = mp_thermo->Y();
        double rho = mp_thermo->density();
        double me = mp_thermo->speciesMw(0);
        double T = mp_thermo->T();
        double Tv = mp_thermo->Tv();

        double tau = compute_tau_ET();
        return 3.E0*RU*rho*p_Y[0]*(T-Tv)/(2.0*me*tau);
    }
    
private:
    const Thermodynamics::Thermodynamics* mp_thermo;
    Transport::Transport* mp_transport;
    Transport::CollisionDB* mp_collisions;

    double const compute_tau_ET();
};

/**
 *  * Represents a coupling between chemistry and electrons.
 *   */
class OmegaCE : public TransferModel
{

public:
    OmegaCE(const Thermodynamics::Thermodynamics& thermo, Kinetics::Kinetics& kinetics){
        mp_thermo = &thermo;
        mp_kinetics = &kinetics;
        m_ns = thermo.nSpecies();
        mp_wrk1 = new double [m_ns];
        mp_wrk2 = new double [m_ns];
    };

    ~OmegaCE()
    {
        delete [] mp_wrk1;
        delete [] mp_wrk2;
    };

    double source()
    {
      static int i_transfer_model = 0;
      switch (i_transfer_model){
          case 0:
              return compute_source_Candler();
              break;
          default:
              std::cerr << "The selected Electron-Chemistry-Electron model is not implemented yet";
              return 0.0;
      }
    };

private:
    const Thermodynamics::Thermodynamics* mp_thermo;
    Kinetics::Kinetics* mp_kinetics;
    int m_ns;
    double* mp_wrk1;
    double* mp_wrk2;

    double const compute_source_Candler();
};


class OmegaEV : public TransferModel
{
  
};

class OmegaI : public TransferModel
{

public:
    OmegaI(const Thermodynamics::Thermodynamics& thermo, Kinetics::Kinetics& kinetics) {
        mp_thermo = &thermo;
        mp_kinetics = &kinetics;
        m_ns = thermo.nSpecies();
        m_nr = kinetics.nReactions();
        mp_hf = new double [m_ns];
        mp_h = new double [m_ns];
        mp_rate = new double [m_nr];
        mp_delta = new double [m_nr];
        for(int i=0; i<m_nr; ++i) {
            if ( (kinetics.reactions()[i].type() == Kinetics::IONIZATION_E)
              || (kinetics.reactions()[i].type() == Kinetics::ION_RECOMBINATION_E)
              || (kinetics.reactions()[i].type() == Kinetics::DISSOCIATION_E)
              || (kinetics.reactions()[i].type() == Kinetics::RECOMBINATION_E))
                 m_rId.push_back(i);
        }
    };

    ~OmegaI() {
        delete [] mp_hf;
        delete [] mp_h;
        delete [] mp_rate;
        delete [] mp_delta;
    };

    /**
      * Computes the Electron-Impact reactions heat Generation in \f$ [J/(m^3\cdot s)] \f$
      * which acts as a shrink to the free electron energy equation.
      *
      * \f[ \Omega^{I} = - \sum_{r \in \mathcal{R}} \Delta h_r \xi_r \f]
      * 
      * \f[ \mathcal{R} \f] denotes the set of electron impact reactions.
      * \f[ \Delta h_r \f] is the reaction enthalpy \f[ [J/mol] \f]
      * \f[ \xi_r \f] is the molar rate of progress \f[ [mol/(m^3\cdot s)] \f]
      *
      */

     double source()
     {
        // Get Formation enthalpy
        mp_thermo->speciesHOverRT(mp_h, NULL, NULL, NULL, NULL, mp_hf);
        for (int i=0; i< m_ns-1; ++i)
            mp_hf[i]*= RU*mp_thermo->T();

        // Get reaction enthapies
        std::fill(mp_delta, mp_delta+m_nr, 0.0);
        mp_kinetics->getReactionDelta(mp_hf,mp_delta);

        // Get molar rates of progress
        mp_kinetics->netRatesOfProgress(mp_rate);

        double src=0;
        for (int i=0; i< m_rId.size(); i++) 
            src -= mp_delta[m_rId[i]]*mp_rate[m_rId[i]];
        return src;
      };

private:
    const Thermodynamics::Thermodynamics* mp_thermo;
    Kinetics::Kinetics* mp_kinetics;
    int m_ns;
    int m_nr;
    std::vector<int> m_rId;
    double* mp_hf;
    double* mp_h;
    double* mp_rate;
    double* mp_delta;
};

/// @}

    } // namespace Transfer
} // namespace Mutation

#endif // TRANSFER_TRANSFER_MODEL_H
