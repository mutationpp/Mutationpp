#ifndef TRANSFER_TRANSFER_MODEL_H
#define TRANSFER_TRANSFER_MODEL_H

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "Transport.h"
#include "MillikanWhite.h"

namespace Mutation {
    namespace Transfer {

class TransferModel
{
public:
    virtual ~TransferModel() { }
    virtual void source(double* const p_source) = 0; 
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

            n_species                 = thermo.nSpecies();
            m_transfer_nHeavy         = thermo.nHeavy();
            m_transfer_offset         = thermo.hasElectrons() ? 1 : 0;
    
	    mp_Mw = new double [n_species];
	    
            for(int i_Mw = 0; i_Mw < n_species; ++i_Mw){
                mp_Mw[i_Mw] = mp_thermo->speciesMw(i_Mw);               
            }
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
     * More information about the above model can be found in: 
     * 
     * C. Park, Review of chemical kinetics problems of future NASA missions, I: Earth entries,
     * Journal of Thermophysics and Heat Transfer, 1993, 7(3):385.
     * 
     * R. N. Schwarz, Z. I. Slawsky, K. F. Herzfeld, Calculation of vibrational relaxation times 
     * in gases, Journal of Chemical Physics, 1952, 20:1591.
     * 
     */
    
    void source(double* const p_source){ 
      
    *p_source = 0.E0;
    const double * p_transfer_Y = mp_thermo->Y();
    double m_transfer_rho = mp_thermo->density();
    // Compute the entalhpy here
    double m_transfer_T = mp_thermo->T();
    double m_transfer_Tv = mp_thermo->Tv();
    
    double h_transfer_total[n_species];  
    double h_transfer_vib_eq[n_species];
    double h_transfer_vib_neq[n_species];
    double mp_work[n_species]; //CHECK IF THIS IS NEEDED
    
    mp_thermo->speciesHOverRT(m_transfer_T, m_transfer_T, m_transfer_T, 
                              m_transfer_T, m_transfer_T, h_transfer_total, mp_work, NULL, 
                              h_transfer_vib_eq, NULL, NULL);

    mp_thermo->speciesHOverRT(m_transfer_T, m_transfer_Tv, m_transfer_T, 
                              m_transfer_Tv, m_transfer_Tv, h_transfer_total, mp_work, NULL, 
                              h_transfer_vib_neq, NULL, NULL);
    
    int i_no_vibrator = 0;
    for (int i_vibrator = 0; i_vibrator-i_no_vibrator < m_mw.nVibrators(); ++i_vibrator){
        if(mp_thermo->species(i_vibrator).type() == Mutation::Thermodynamics::ATOM){
            i_no_vibrator++; 
	} 
        else {
            *p_source += p_transfer_Y[i_vibrator] * m_transfer_rho *RU * m_transfer_T /mp_Mw[i_vibrator] * (h_transfer_vib_eq[i_vibrator] - h_transfer_vib_neq[i_vibrator]) /compute_tau_VT_m(i_vibrator-i_no_vibrator);
	}
    }
	
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
    int n_species;
    int m_transfer_nHeavy; 
    int m_transfer_offset;
    
    double* mp_Mw;
    
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
        p_thermo = &thermo; 
        p_kinetics = &kinetics;
    };
      
    void source(double* const p_source){
      static int i_transfer_model = 0;
      switch (i_transfer_model){
          case 0:
              *p_source = compute_source_Candler();
	      break;
          default:
              std::cerr << "The selected Chemistry-Vibration-Chemistry model is not implemented yet";
      }
    };
    
private:
    const Thermodynamics::Thermodynamics* p_thermo;
    Kinetics::Kinetics* p_kinetics;

    double const compute_source_Candler();

};

class OmegaET : public TransferModel
{
public:
    OmegaET(const Thermodynamics::Thermodynamics& thermo, Transport::Transport& transport){
      
    };
    
    void source(double* const p_source){
        //*p_source = 3.E0/2.E0*KB;
      
    }
    
private:
    double const compute_tau_ET();
};

class OmegaEV : public TransferModel
{
  
};

class OmegaI : public TransferModel
{
  
};

    } // namespace Transfer
} // namespace Mutation

#endif // TRANSFER_TRANSFER_MODEL_H
