#include "StateModel.h"
#include "Utilities.h"
#include <algorithm>

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================


//class EquilTPStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    EquilTPStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//
//    /**
//     * Sets the mixture state by computing the equilibrium composition at the
//     * given temperature and pressure using the default elemental composition.
//     *
//     * @param p_T - pointer to temperature value (K)
//     * @param p_P - pointer to pressure value (Pa)
//     */
//    virtual void setState(const double* const p_T, const double* const p_P)
//    {
//        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_T[0];
//        m_P = p_P[0];
//        m_thermo.equilibriumComposition(m_T, m_P, mp_X);
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    EquilTPStateModel, StateModel> equil_tp("EquilTP");
//
////==============================================================================
//
//class EquilTPXeStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    EquilTPXeStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//
//    /**
//     * Sets the mixture state by computing the equilibrium composition at the
//     * given temperature and pressure using the default elemental composition.
//     *
//     * @param p_T - pointer to temperature value (K)
//     * @param p_P - pointer to pressure value (Pa)
//     */
//    virtual void setState(const double* const p_TP, const double* const p_Xe)
//    {
//        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_TP[0];
//        m_P = p_TP[1];
//        m_thermo.equilibriumComposition(m_T, m_P, p_Xe, mp_X);
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    EquilTPXeStateModel, StateModel> equil_tpxe("EquilTPXe");
//
////==============================================================================
//
//class TPXStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    TPXStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//    
//    /**
//     * Sets the mixture state from mixture temperature, pressure, and species
//     * mole fractions.
//     *
//     * @param p_TP - pointer to {T (K), P (Pa)} array
//     * @param p_X  - pointer to species mole fraction array
//     */
//    virtual void setState(const double* const p_TP, const double* const p_X)
//    {
//        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_TP[0];
//        m_P = p_TP[1];
//        std::copy(p_X, p_X+m_thermo.nSpecies(), mp_X);
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    TPXStateModel, StateModel> tpx("TPX");
//
////==============================================================================
//
//class TPYStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    TPYStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//    
//    /**
//     * Sets the mixture state from mixture temperature, pressure, and species
//     * mass fractions.
//     *
//     * @param p_TP - pointer to {T (K), P (Pa)} array
//     * @param p_Y  - pointer to species mass fraction array
//     */
//    virtual void setState(const double* const p_TP, const double* const p_Y)
//    {
//        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_TP[0];
//        m_P = p_TP[1];
//        m_thermo.convert<Y_TO_X>(p_Y, mp_X);
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    TPYStateModel, StateModel> tpy("TPY");
//
////==============================================================================
//
//class TRhoiStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    TRhoiStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//    
//    /**
//     * Sets the mixture state from mixture temperature and species densities.
//     *
//     * @param p_T    - pointer to mixture temperature
//     * @param p_rhoi - pointer to species densities array
//     */
//    virtual void setState(const double* const p_T, const double* const p_rhoi)
//    {
////        cout << "T,Rhoi : setState(" << p_T[0];
////        for (int i = 0; i < m_thermo.nSpecies(); ++i)
////            cout << ", " << p_rhoi[i];
////        cout << ")" << endl;
//        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_T[0];
//        
//        m_P = 0.0;
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            m_P += p_rhoi[i];
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            mp_X[i] = p_rhoi[i] / m_P;
//        m_P = m_thermo.pressure(m_T, m_P, mp_X);
//        
//        m_thermo.convert<Y_TO_X>(mp_X, mp_X);
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    TRhoiStateModel, StateModel> trhoi("TRhoi");
//
////==============================================================================
//
//class TTvRhoiStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    TTvRhoiStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//    
//    /**
//     * Sets the mixture state from mixture temperature, pressure, and species
//     * mass fractions.
//     *
//     * @param p_T    - pointer to {T (K), Tv (K)} array
//     * @param p_rhoi - pointer to species densities array
//     */
//    virtual void setState(const double* const p_T, const double* const p_rhoi)
//    {
//        m_T = m_Tr = p_T[0];
//        m_Tv = m_Tel = m_Te = p_T[1];
//        
//        m_P = 0.0;
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            m_P += p_rhoi[i];
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            mp_X[i] = p_rhoi[i] / m_P;
//        m_P = m_thermo.pressure(m_T, m_P, mp_X);
//        
//        m_thermo.convert<Y_TO_X>(mp_X, mp_X);
//    }
//    
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    TTvRhoiStateModel, StateModel> ttvrhoi("TTvRhoi");
//
////==============================================================================
//
//class TTvTeRhoiStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    TTvTeRhoiStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//    
//    /**
//     * Sets the mixture state from mixture temperature, pressure, and species
//     * mass fractions.
//     *
//     * @param p_T    - pointer to {T (K), Tv (K), Te (K)} array
//     * @param p_rhoi - pointer to species densities array
//     */
//    virtual void setState(const double* const p_T, const double* const p_rhoi)
//    {
//        m_T = m_Tr = p_T[0];
//        m_Tv = p_T[1];
//        m_Tel = m_Te = p_T[2];
//        
//        m_P = 0.0;
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            m_P += p_rhoi[i];
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            mp_X[i] = p_rhoi[i] / m_P;
//        m_P = m_thermo.pressure(m_T, m_P, mp_X);
//        
//        m_thermo.convert<Y_TO_X>(mp_X, mp_X);
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    TTvTeRhoiStateModel, StateModel> ttvterhoi("TTvTeRhoi");
//
////==============================================================================
//// This is the state model for the stagnation line code single temperature
//
//class RhoiRhoEStateModel : public StateModel
//{
//public:
//    /**
//     * Constructor.
//     */
//    RhoiRhoEStateModel(const Thermodynamics& thermo)
//        : StateModel(thermo)
//    { }
//    
//    /**
//     * Sets the mixture state from mixture partial densities and internal energy
//     *
//     * @param p_rhoi - pointer to species densities array
//     * @param p_E - Internal energy
//     */
//    virtual void setState(const double* const p_rhoi, const double* const p_rhoE)
//    {
//    int nb_temp = 1;
//    int nb_ns = m_thermo.nSpecies(); 
//    // Getting temperature from p_E. Newton iteration method. It should be implemented somewhere else
//    //STEPS: 
//    //0) Initial guess of temperature
//    // Initialize more things in general
//    double T = 3000.e0;
//    double dT = 0;
//    double residual = 1.e0; //Initialize
//    double f, fp;
//    
//    double * p_h = new double[nb_ns];
//    double * p_cp = new double[nb_ns];
//    double * p_Ri = new double[nb_ns];
//
//   while(residual > 1.e-8){
//
//     //2) get hi's for ei mix.speciesHOverRT(double T, double* const h)
//
//     m_thermo.speciesHOverRT(T, p_h);
//     m_thermo.speciesCpOverR(T, p_cp);
// 
//     //3) get Ri
//     //4) get f, fp
//     f = 0; fp = 0;
//
//     for(int ic_1=0; ic_1 < nb_ns; ic_1++)
//     {
//       p_Ri[ic_1] = (RU/m_thermo.speciesMw(ic_1)); //check
//       p_h[ic_1] *= T * (p_Ri[ic_1]); //check
//       p_cp[ic_1] *= p_Ri[ic_1]; //check
//       fp += p_rhoi[ic_1] * (p_cp[ic_1] - p_Ri[ic_1]);
//       f += p_rhoi[ic_1] * (p_h[ic_1] - T * p_Ri[ic_1]);
//     }
//
//     //5) dT = (f - e)/fp
//     //6) T = T - dT
//     dT = (f - *p_rhoE)/fp;
//     T = T - dT;
//
//     //7) compute residual
//     residual = std::abs(dT/T);
//
//     m_T = m_Tr = m_Tv = m_Tel = m_Te = T;
//
//    }
//
//     // COPY-PASTED
//        m_P = 0.0;
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            m_P += p_rhoi[i];
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            mp_X[i] = p_rhoi[i] / m_P;
//        m_P = m_thermo.pressure(m_T, m_P, mp_X);
//        
//        m_thermo.convert<Y_TO_X>(mp_X, mp_X);
//
//
//    delete [] p_h;
//    delete [] p_cp;
//    delete [] p_Ri;
//
//    }
//};
//
//// Register the state model
//Utilities::Config::ObjectProvider<
//    RhoiRhoEStateModel, StateModel> rhoirhoe("RhoiRhoE");
//
//==============================================================================
//============    HERE IS THE NEW CATEGORY OF STATE MODELS    ==================
//==============================================================================

// This is the state model for the stagnation line code single temperature

class NonequilibriumSingleT : public StateModel
{
//Determines the state variables  in order to characterize the mixture via set_State function
// state_variable = 0 --> Conservative variables; partial mass densities and total energy density
// state_variable = 1 --> Partial densities and temperature
int state_variable;

public:
    /**
     * Constructor.
     */
    NonequilibriumSingleT(const Thermodynamics& thermo)
        : StateModel(thermo)
    { 
// Default state_variable initialized
    state_variable = 0;
    }
    
    virtual void setState(const double* const p_mass, const double* const p_energy, int i_state)
    {
      state_variable = i_state;

      switch (state_variable)
      {
        case 0:
        {
          // Default model --> rhoie, rhoE
        
          int nb_temp = nEnergyEq();
          int nb_ns = m_thermo.nSpecies(); 

          double T = 1000.e0;
          double dT = 0;
          double residual = 1.e0;
          double f, fp;
    
          double * p_h = new double[nb_ns];
          double * p_cp = new double[nb_ns];
          double * p_Ri = new double[nb_ns];

          while(residual > 1.e-8)
          {
            m_thermo.speciesHOverRT(T, p_h);
            m_thermo.speciesCpOverR(T, p_cp);
            f = 0; fp = 0;

            for(int ic_1=0; ic_1 < nb_ns; ic_1++)
            {
              p_Ri[ic_1] = (RU/m_thermo.speciesMw(ic_1));
              p_h[ic_1] *= T * (p_Ri[ic_1]);
              p_cp[ic_1] *= p_Ri[ic_1];
              fp += p_mass[ic_1] * (p_cp[ic_1] - p_Ri[ic_1]);
              f += p_mass[ic_1] * (p_h[ic_1] - T * p_Ri[ic_1]);
            }

            dT = (f - *p_energy)/fp;
            T = T - dT;

            residual = std::abs(dT/T);

            m_T = m_Tr = m_Tv = m_Tel = m_Te = T;

          }

          m_P = 0.0;
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_mass[i];
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_mass[i] / m_P;
          m_P = m_thermo.pressure(m_T, m_P, mp_X);

          m_thermo.convert<Y_TO_X>(mp_X, mp_X);

          delete [] p_h;
          delete [] p_cp;
          delete [] p_Ri;

          break;
        }
        case 1:
        {
          // rhoi, T

          m_T = m_Tr = m_Tv = m_Tel = m_Te = *p_energy;
        
          m_P = 0.0;
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_mass[i];
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_mass[i] / m_P;
          m_P = m_thermo.pressure(m_T, m_P, mp_X);
          m_thermo.convert<Y_TO_X>(mp_X, mp_X);

          break;
        }
        default:
        {
          std::cout << "The setState model requested has not been implemented yet." << endl;
        }

      }

    }
    
    virtual int nEnergyEq()
    {
      return(1);
    }

    virtual int nVibEq()
    {
      return(0);
    }

    virtual int nRotEq()
    {
      return(0);
    }

    virtual int nElEq()
    {
      return(0);
    }

};

// Register the state model
Utilities::Config::ObjectProvider<
    NonequilibriumSingleT, StateModel> nonequilibriumsinglet("NonequilibriumSingleT");

// This is the state model for the stagnation line code two temperature

class NonequilibriumTTv : public StateModel
{
//Determines the state variables  in order to characterize the mixture via set_State function
// state_variable = 0 --> Conservative variables; partial mass densities, total energy density and vib energy
// state_variable = 1 --> Partial densities and temperatures T and Tv
int state_variable;

public:
    /**
     * Constructor.
     */
    NonequilibriumTTv(const Thermodynamics& thermo)
        : StateModel(thermo)
    { 
// Default state_variable initialized
    state_variable = 0;
    }
    
    virtual void setState(const double* const p_mass, const double* const p_energy, int i_state)
    {
      state_variable = i_state;

      switch (state_variable)
      {
        case 0:
        {
          // Default model --> rhoie, rhoE

          int nb_temp = nEnergyEq();
          int nb_ns = m_thermo.nSpecies(); 
 
          double T=300.e0;
          double dT = 0.e0;
          double Tv=300.e0;
          double Tv_trial = 0;
          double dTv = 0;
          double residual = 1.e0;
          double residual_temp = 1.e0;
          double cptr, hf, fv, fpv;
          double f, fp;

          bool b_in;

          double * p_h = new double[nb_ns];
          double * p_ht = new double[nb_ns];
          double * p_hr = new double[nb_ns];
          double * p_hv = new double[nb_ns];
          double * p_he = new double[nb_ns];
          double * p_hf = new double[nb_ns];
          double * p_cp = new double[nb_ns];
          double * p_cpt = new double[nb_ns];
          double * p_cpr = new double[nb_ns];
          double * p_cpv = new double[nb_ns];
          double * p_cpe = new double[nb_ns];
          double * p_Ri = new double[nb_ns];

          for(int ic_1=0; ic_1 < nb_ns; ++ic_1){p_Ri[ic_1] = (RU/m_thermo.speciesMw(ic_1));}

          // Computing Tv

            m_T = m_Tr = T;
            m_Tv = m_Tel = m_Te = Tv;

          while(residual > 1.e-8){
            m_thermo.speciesHOverRT(p_h, p_ht, p_hr, p_hv, p_he, p_hf);
            m_thermo.speciesCpOverR(T, Tv, T, Tv, Tv, p_cp, p_cpt, p_cpr, p_cpv, p_cpe);

            fv = 0; fpv = 0;

            for(int ic_1=0; ic_1 < nb_ns; ic_1++)
            {
              p_hv[ic_1] *= T * (p_Ri[ic_1]);
              p_cpv[ic_1] *= p_Ri[ic_1];

              fpv += p_mass[ic_1] * (p_cpv[ic_1]);
              fv += p_mass[ic_1] * (p_hv[ic_1]);
            }

            dTv = (fv - p_energy[1])/fpv;
            Tv_trial = Tv - dTv;
            residual = std::abs(fv - p_energy[1]);

            b_in = true;

//          Test for residual
            double lambda = 1.0;             
            while(b_in){
              m_Tv = m_Tel = m_Te = Tv_trial;

              m_thermo.speciesHOverRT(p_h, p_ht, p_hr, p_hv, p_he, p_hf);

              fv =0;

              for(int ic_1=0; ic_1 < nb_ns; ic_1++)
              {
                p_hv[ic_1] *= T * (p_Ri[ic_1]);

                fv += p_mass[ic_1] * (p_hv[ic_1]);
              }

              residual_temp = std::abs(fv - p_energy[1]);

              if(residual_temp > (1.0-1.0e-4*lambda)*residual){
                  lambda *= 0.5;
                  Tv_trial = Tv - lambda*dTv;
              }
              else{
                Tv = Tv_trial;
                b_in = false;
              }
            }
            m_Tv = m_Tel = m_Te = Tv;
          }

          // Imposing equilibrium as initial guess

          m_T = m_Tr = T = Tv;

          // Computing T

          while(residual > 1.e-8)
          {
            m_thermo.speciesHOverRT(T, p_h);
            m_thermo.speciesCpOverR(T, p_cp);
            f = 0; fp = 0;

            for(int ic_1=0; ic_1 < nb_ns; ic_1++)
            {
              p_h[ic_1] *= T * (p_Ri[ic_1]);
              p_cp[ic_1] *= p_Ri[ic_1];
              fp += p_mass[ic_1] * (p_cp[ic_1] - p_Ri[ic_1]);
              f += p_mass[ic_1] * (p_h[ic_1] - T * p_Ri[ic_1]);
            }

            dT = (f - *p_energy)/fp;
            T = T - dT;

            // Bad residual fix me

            residual = std::abs(dT/T);

          }

//          cout << T << " " << Tv << endl;

          // Setting the final state

          m_T = m_Tr = T;
          m_Tv = m_Tel = m_Te = Tv;

          m_P = 0.0;
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_mass[i];
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_mass[i] / m_P;
          m_P = m_thermo.pressure(m_T, m_P, mp_X);

          m_thermo.convert<Y_TO_X>(mp_X, mp_X);

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//          int nb_temp = nEnergyEq();
//          int nb_ns = m_thermo.nSpecies(); 
//
//          double T; //No need to be initialized
//          double Tv = 1000.e0;
//          double Tv_trial = 0;
//          double dTv = 0;
//          double residual = 1.e0;
//          double residual_temp = 1.e0;
//          double cptr, hf, fv, fpv;
//
//          bool b_in;
//    
//          double * p_h = new double[nb_ns];
//          double * p_ht = new double[nb_ns];
//          double * p_hr = new double[nb_ns];
//          double * p_hv = new double[nb_ns];
//          double * p_he = new double[nb_ns];
//          double * p_hf = new double[nb_ns];
//          double * p_cp = new double[nb_ns];
//          double * p_cpt = new double[nb_ns];
//          double * p_cpr = new double[nb_ns];
//          double * p_cpv = new double[nb_ns];
//          double * p_cpe = new double[nb_ns];
//          double * p_Ri = new double[nb_ns];
//
//          cptr = 0;
//
//          for(int ic_1=0; ic_1 < nb_ns; ++ic_1){p_Ri[ic_1] = (RU/m_thermo.speciesMw(ic_1));}
//
//          m_T = m_Tr = m_Tv = m_Tel = m_Te = Tv;
//
//          m_thermo.speciesHOverRT(p_h, NULL, NULL, NULL, NULL, p_hf);
//
//          for(int ic_1=0; ic_1 < nb_ns; ++ic_1){
//          hf += p_mass[ic_1] * p_hf[ic_1] * Tv * p_Ri[ic_1];
//          }
//
////        Translational Temperature
//          for(int ic_1=0; ic_1 < nb_ns; ++ic_1){
//            cptr += p_mass[ic_1] * 1.5e0 * p_Ri[ic_1];
//            if(m_thermo.species(ic_1).type() == MOLECULE ? 1 : 0) {cptr += p_mass[ic_1] * p_Ri[ic_1];} //Only for linear molecules
//            }
//
//          T = (p_energy[0] - p_energy[1] - hf)/(cptr);
//          Tv = T;
//
//          residual = 1.e0;
//
//          m_T = m_Tr = T;
//          m_Tv = m_Tel = m_Te = Tv;
//
////        Newton loop for Tv 
//
//          while(residual > 1.e-8){
//            m_thermo.speciesHOverRT(p_h, p_ht, p_hr, p_hv, p_he, p_hf);
//            m_thermo.speciesCpOverR(T, Tv, T, Tv, Tv, p_cp, p_cpt, p_cpr, p_cpv, p_cpe);
//
//            fv = 0; fpv = 0;
//
////          SUPER SOS: MAKE IT GENERAL SO THAT IT CAN TAKE INTO ACCOUNT FREE ELECTRONS ALSO
//
//            for(int ic_1=0; ic_1 < nb_ns; ic_1++)
//            {
//              p_hv[ic_1] *= T * (p_Ri[ic_1]);
//              p_he[ic_1] *= T * (p_Ri[ic_1]);
//              p_cpv[ic_1] *= p_Ri[ic_1];
//              p_cpe[ic_1] *= p_Ri[ic_1];
//
//              fpv += p_mass[ic_1] * (p_cpv[ic_1] + p_cpe[ic_1]);
//              fv += p_mass[ic_1] * (p_hv[ic_1] + p_he[ic_1] );
//            }
//
//            dTv = (fv - p_energy[1])/fpv;
//            Tv_trial = Tv - dTv;
//            residual = std::abs(fv - p_energy[1]);
//
//            b_in = true;
//
////          Test for residual
//            double lambda = 1.0;             
//            while(b_in){
//              m_T = m_Tr = T;
//              m_Tv = m_Tel = m_Te = Tv_trial;
//
//              m_thermo.speciesHOverRT(p_h, p_ht, p_hr, p_hv, p_he, p_hf);
//
//              fv =0;
//
//              for(int ic_1=0; ic_1 < nb_ns; ic_1++)
//              {
//                p_hv[ic_1] *= T * (p_Ri[ic_1]);
//                p_he[ic_1] *= T * (p_Ri[ic_1]);
//
//                fv += p_mass[ic_1] * (p_hv[ic_1] + p_he[ic_1]);
//              }
//
//              residual_temp = std::abs(fv - p_energy[1]);
//
//              if(residual_temp > (1.0-1.0e-4*lambda)*residual){
//                  lambda *= 0.5;
//                  Tv_trial = Tv - lambda*dTv;
//              }
//              else{
//                Tv = Tv_trial;
//                b_in = false;
//              }
//            }
//            m_Tv = m_Tel = m_Te = Tv;
//          }
//
//          m_P = 0.0;
//          for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            m_P += p_mass[i];
//          for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            mp_X[i] = p_mass[i] / m_P;
//          m_P = m_thermo.pressure(m_T, m_P, mp_X);
//
//          m_thermo.convert<Y_TO_X>(mp_X, mp_X);
//
          delete [] p_h;
          delete [] p_ht;
          delete [] p_hr;
          delete [] p_hv;
          delete [] p_he;
          delete [] p_hf;
          delete [] p_cp;
          delete [] p_cpt;
          delete [] p_cpr;
          delete [] p_cpv;
          delete [] p_cpe;
          delete [] p_Ri;

          break;
        }
        case 1:
        {
          // rhoi, T

          m_T = m_Tr = p_energy[0];
          m_Tv = m_Tel = m_Te = p_energy[1];
        
          m_P = 0.0;
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_mass[i];
          for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_mass[i] / m_P;
          m_P = m_thermo.pressure(m_T, m_P, mp_X);
          m_thermo.convert<Y_TO_X>(mp_X, mp_X);

          break;
        }
        default:
        {
          std::cout << "The setState model requested has not been implemented yet." << endl;
        }

      }

    }
    
    virtual int nEnergyEq()
    {
      return(2);
    }

    virtual int nVibEq()
    {
      return(1);
    }

    virtual int nRotEq()
    {
      return(0);
    }

    virtual int nElEq()
    {
      return(0);
    }

};

// Register the state model
Utilities::Config::ObjectProvider<
    NonequilibriumTTv, StateModel> nonequilibriumttv("NonequilibriumTTv");


    } // namespace Thermodynamics
} // namespace Mutation


