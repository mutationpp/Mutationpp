#include "StateModel.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * Chemical Non-equilibrium 2-T StateModel. Thermal nonequilibrium with
 * T = Tr and Tv = Tel = Te with finite rate chemistry.
 */
class ChemNonEqTTvStateModel : public StateModel
{
public:
  
     //Mutation::Transfer::TransferModel* p_transfer_model;

    ChemNonEqTTvStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 2, thermo.nSpecies())
    {
        mp_work1 = new double [thermo.nSpecies()];
        mp_work2 = new double [thermo.nSpecies()];
        mp_work3 = new double [thermo.nSpecies()];
        mp_work4 = new double [thermo.nSpecies()];
    }

    ~ChemNonEqTTvStateModel()
    {
        delete [] mp_work1;
        delete [] mp_work2;
        delete [] mp_work3;
        delete [] mp_work4;
        delete mp_omega_VT;
        delete mp_omega_CV;
        delete mp_omega_CE;
//        delete mp_omega_EV;
//        delete mp_omega_ET;
//        delete mp_omega_I;
        //delete p_transfer_model;
    }

    /**
     * Sets the current state of the mixture.  Variable sets can be the
     * following:
     *   0: {species densities}, {total energy density, vib. energy density}
     *   1: {species densities} and {T, Tv}
     */
    void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0)
    {
        const int ns = m_thermo.nSpecies();

        // Compute the species concentrations which are used through out this
        // method regardless of variable set 
        double conc = 0.0;
        for (int i = 0; i < ns; ++i){
            // Check that species densities are at least positive
            assert(p_mass[i] >= 0.0);
            mp_X[i] = p_mass[i] / m_thermo.speciesMw(i);
            conc += mp_X[i];
        }

        // Compute the temperatures and make sure the variable set is implemented
        switch (vars) {
        case 0: {
            // First step is to solve one nonlinear equation to get Tv from the
            // vibrational energy density equation

            getTFromRhoE(
                Cpv(m_thermo), Hv(m_thermo), p_energy[1], m_Tv, mp_work1, 0.0);

            // Next compute the temperature using the total energy density

            m_Tel = m_Te = m_Tv;

            getTFromRhoE(
                Cp(m_thermo), H(m_thermo, m_Tv), p_energy[0], m_T, mp_work1, -conc);

            break;
        }
        case 1:
            // Check that temperature is at least positive
            assert(p_energy[0] > 0.0);
            assert(p_energy[1] > 0.0);
            m_T = p_energy[0];
            m_Tv = p_energy[1];
            break;
        default:
            cout << "Variable-set " << vars << " not implemented in StateModel!" << endl;
            exit(1);
        }

        // Set Tr, Tel, and Te
        m_Tr = m_T;
        m_Tel = m_Te = m_Tv;

        // Compute the pressure and species mole fractions from T and rho_i
        for (int i = 0; i < ns; ++i)        
            mp_X[i] /= conc;
        m_P = RU * m_T * conc;

        /* REMOVE ME
        p_transfer_model = new Mutation::Transfer::TransferModelVT(m_thermo); // I include a copy of thermo
        double output;
        p_transfer_model->source(&output);

        // std::cout << "Source VT = " << output << endl;
        */

    }
    
    void initializeTransferModel(
        Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        Mutation::Kinetics::Kinetics& kinetics){
      
       mp_omega_VT = new Mutation::Transfer::OmegaVT(thermo);
       mp_omega_CV = new Mutation::Transfer::OmegaCV(thermo, kinetics);
       mp_omega_CE = new Mutation::Transfer::OmegaCE(thermo, kinetics);
       //mp_omega_ET = new Mutation::Transfer::OmegaET();
       //mp_omega_EV = new Mutation::Transfer::OmegaEV();
       //mp_omega_I  = new Mutation::Transfer::OmegaI();
       
    }
    
    void energyTransferSource(double* const omega)
    {
            double m_work;
      
            omega[0] = 0.E0;
            mp_omega_VT->source(&m_work); 
            omega[0] += m_work;
            mp_omega_CV->source(&m_work); 
            omega[0] += m_work;
            mp_omega_CE->source(&m_work); 
            omega[0] += m_work;
//            mp_omega_ET->source(&m_work); 
//            omega[0] += m_work;
//            mp_omega_EV->source(&m_work); 
//            omega[0] += m_work;
//            mp_omega_I->source(&m_work); 
//            omega[0] += m_work;
    }
    
    void getTemperatures(double* const p_T) const {
        p_T[0] = T();
        p_T[1] = Tv();
    }
    
    void getEnergyMass(double* const p_e) {
    
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work1, NULL, NULL, mp_work2, mp_work3, NULL);
        double species_Mw[n_species];
        
        for(int i_energy_mass = 0; i_energy_mass < m_thermo.nSpecies(); ++i_energy_mass){
            species_Mw[i_energy_mass] = m_thermo.speciesMw(i_energy_mass);
            p_e[i_energy_mass] = mp_work1[i_energy_mass] * m_T * RU /species_Mw[i_energy_mass] - RU/species_Mw[i_energy_mass] *m_T;
        }
        for(int i_energy_mass = 0; i_energy_mass < m_thermo.nSpecies(); ++i_energy_mass){
            p_e[i_energy_mass+n_species] = (mp_work2[i_energy_mass] + mp_work3[i_energy_mass]) 
                                            * m_T * RU / species_Mw[i_energy_mass];
        }
        
    }

    void getEnthalpyMass(double* const p_h) {
        
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work1, NULL, NULL, mp_work2, mp_work3, NULL);
        double species_Mw[n_species];

        for(int i_enthalpy_mass = 0; i_enthalpy_mass < n_species; ++i_enthalpy_mass){
            species_Mw[i_enthalpy_mass] = m_thermo.speciesMw(i_enthalpy_mass);
            p_h[i_enthalpy_mass] = mp_work1[i_enthalpy_mass] * m_T * RU /species_Mw[i_enthalpy_mass];
        }
        for(int i_enthalpy_mass = 0; i_enthalpy_mass < n_species; ++i_enthalpy_mass){
            p_h[i_enthalpy_mass+n_species] = (mp_work2[i_enthalpy_mass] + mp_work3[i_enthalpy_mass])
                                              * m_T * RU / species_Mw[i_enthalpy_mass];
        }

    }
    
    void getCpMass(double* const p_Cp) {
        
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesCpOverR(m_T, m_Te, m_Tr, m_Tv, m_Tel, NULL, mp_work1, mp_work2, mp_work3, mp_work4);
        double species_Mw[n_species];

        for(int i_getCp = 0; i_getCp < n_species; ++i_getCp){
            species_Mw[i_getCp] = m_thermo.speciesMw(i_getCp);
            p_Cp[i_getCp] = (mp_work1[i_getCp] + mp_work2[i_getCp]) *RU/species_Mw[i_getCp];
        }
        for(int i_getCp = 0; i_getCp < n_species; ++i_getCp){
            p_Cp[i_getCp + n_species] = (mp_work3[i_getCp] + mp_work4[i_getCp])*RU/species_Mw[i_getCp];
        }
        
    }

    void getTagModes(int* const p_tag) {

      p_tag[0] = 1; p_tag[5] = 0; // Heavy translation
      p_tag[1] = 0; p_tag[6] = 1; // Electron translation
      p_tag[2] = 1; p_tag[7] = 0; // Rotation excitation
      p_tag[3] = 0; p_tag[8] = 1; // Vibration excitation
      p_tag[4] = 0; p_tag[9] = 1; // Electronic excitation

    }
    
    private:

    /**
     * Small helper class which provides wrapper to
     * Thermodynamics::speciesCpOverR() to be used with getTFromRhoE().
     */
    
    class Cp {
    public:
        Cp(const Thermodynamics& t) : thermo(t) {
             p_cpt = new double [thermo.nSpecies()]; 
             p_cpr = new double [thermo.nSpecies()]; 
             p_cpv = new double [thermo.nSpecies()];
             p_cpel = new double [thermo.nSpecies()];
        }
        void operator () (double T, double* const cp) const {
          thermo.speciesCpOverR(T, T, T, T, T, NULL, p_cpt, p_cpr, NULL, NULL);
          int n_species = thermo.nSpecies();
          for(int iCp = 0; iCp < n_species; ++iCp) {
            cp[iCp] = p_cpt[iCp]+p_cpr[iCp];
          }
        }

    private:
        const Thermodynamics& thermo;
	double* p_cpt;
	double* p_cpr;
	double* p_cpv;
	double* p_cpel;
    }; // class Cp
    
    class Cpv {
    public:
        Cpv(const Thermodynamics& t) : thermo(t) {
             p_cpv = new double [thermo.nSpecies()];
             p_cpel = new double [thermo.nSpecies()];
        }
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(T, T, T, T, T, NULL, NULL, NULL, p_cpv, p_cpel);
            for(int iCpv = 0; iCpv < thermo.nSpecies(); ++iCpv){
                cp[iCpv] = p_cpv[iCpv] + p_cpel[iCpv];
            }
        }
    private:
        const Thermodynamics& thermo;
	double* p_cpv;
	double* p_cpel;
    }; // class Cpv

    /**
     * Small helper classes which provide wrapper to
     * Thermodynamics::speciesHOverRT() to be used with getTFromRhoE().
     */
    
    class H {
    public:
        H(const Thermodynamics& t, const double& T_vibration) : thermo(t) {
            Tv = T_vibration;
        }
        void operator () (double T, double* const h) const {
          thermo.speciesHOverRT(T, Tv, T, Tv, Tv, h, NULL, NULL, NULL, NULL, NULL);
        }
    private:
        const Thermodynamics& thermo;
        double Tv;
    }; // class H
    
    class Hv {
    public:
        Hv(const Thermodynamics& t) : thermo(t) { 
             p_hv  = new double [thermo.nSpecies()]; 
             p_hel = new double [thermo.nSpecies()]; 
        }
        void operator () (double T, double* const h) const {
            thermo.speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, p_hv, p_hel, NULL);
            for(int ihv = 0; ihv < thermo.nSpecies(); ++ihv){
                h[ihv] = p_hv[ihv] + p_hel[ihv];
            }
        }
    private:
        const Thermodynamics& thermo;
        double* p_hv;
        double* p_hel;
    }; // class Hv

private:
    double* mp_work1;
    double* mp_work2;
    double* mp_work3;
    double* mp_work4;
    Mutation::Transfer::TransferModel* mp_omega_VT;
    Mutation::Transfer::TransferModel* mp_omega_CV;
    Mutation::Transfer::TransferModel* mp_omega_CE;
//    Mutation::Transfer::TransferModel* mp_omega_EV;
//    Mutation::Transfer::TransferModel* mp_omega_ET;
//    Mutation::Transfer::TransferModel* mp_omega_I;
    
}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqTTvStateModel, StateModel> chem_non_eq_TTv("ChemNonEqTTv");

    } // namespace Thermodynamics
} // namespace Mutation
