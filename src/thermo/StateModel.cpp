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
//    // Get it robustly //!!!//
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
    
    virtual void setState(const double* const p_mass, const double* const p_energy, int state)
    {
      state_variable = state;

      switch (state_variable)
      {
      case 0:
      {
        // Default model --> rhoie, rhoE
        
        // Get it robustly //!!!//
        int nb_temp = 1;
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
//
        break;
      }
      default:
      {
        std::cout << "The setState model requested has not been implemented yet." << endl;
      }

    }

    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    NonequilibriumSingleT, StateModel> nonequilibriumsinglet("NonequilibriumSingleT");


    } // namespace Thermodynamics
} // namespace Mutation


