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
        mp_work = new double [thermo.nSpecies()];
    }

    ~ChemNonEqTTvStateModel()
    {
        delete [] mp_work;
        delete mp_omega_VT;
        delete mp_omega_CV;
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

//            m_T = m_Tr = 300.e0;
//            m_Tv = m_Te = m_Tv = 300.e0;

            getTFromRhoE(
                Cpv(m_thermo), Hv(m_thermo), p_energy[1], m_Tv, mp_work, 0.0);

            // Next compute the temperature using the total energy density

            m_Tel = m_Te = m_Tv;

            getTFromRhoE(
                Cp(m_thermo, m_Tv), H(m_thermo, m_Tv), p_energy[0], m_T, mp_work, -conc);

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
//       mp_omega_CV = new Mutation::Transfer::OmegaCV(thermo, kinetics);
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
//            mp_omega_CV->source(&m_work); 
//            omega[0] += m_work;
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
    
    void getEnergyDensities(double* const p_rhoe) {
        std::cerr << "getEnergyDensities()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }
    
    void getEnergyMass(double* const p_e) {
    
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(p_e, mp_work, NULL, p_e+n_species, NULL, NULL);
        double species_Mw[n_species];
        
        for(int i_energy_mass = 0; i_energy_mass < m_thermo.nSpecies(); ++i_energy_mass){
            species_Mw[i_energy_mass] = m_thermo.speciesMw(i_energy_mass);
            p_e[i_energy_mass] = p_e[i_energy_mass] * m_T * RU /species_Mw[i_energy_mass] - RU/species_Mw[i_energy_mass] *m_T;
        }
        for(int i_energy_mass = 0; i_energy_mass < m_thermo.nSpecies(); ++i_energy_mass){
            p_e[i_energy_mass+n_species] *=  m_T * RU / species_Mw[i_energy_mass];
        }
        
    }

    void getEnthalpyDensities(double* const p_rhoh) {
        std::cerr << "getEnthalpyDensities()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }


    void getEnthalpyMass(double* const p_h) {
        
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(p_h, mp_work, NULL, p_h+2, NULL, NULL);
        double species_Mw[n_species];

        for(int i_enthalpy_mass = 0; i_enthalpy_mass < n_species; ++i_enthalpy_mass){
            species_Mw[i_enthalpy_mass] = m_thermo.speciesMw(i_enthalpy_mass);
            p_h[i_enthalpy_mass] *= m_T * RU /species_Mw[i_enthalpy_mass];
        }
        for(int i_enthalpy_mass = 0; i_enthalpy_mass < n_species; ++i_enthalpy_mass){
            p_h[i_enthalpy_mass+n_species] *= m_T * RU / species_Mw[i_enthalpy_mass];
        }

    }
    
    void getCp(double* const p_Cp) {
        
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesCpOverR(m_T, m_Te, m_Tr, m_Tv, m_Tel, p_Cp, NULL, NULL, p_Cp+n_species, NULL);
        double species_Mw[n_species];

        for(int i_getCp = 0; i_getCp < n_species; ++i_getCp){
            species_Mw[i_getCp] = m_thermo.speciesMw(i_getCp);
            p_Cp[i_getCp] *= RU/species_Mw[i_getCp];
        }
        for(int i_getCp = 0; i_getCp < n_species; ++i_getCp){
            p_Cp[i_getCp + n_species] *= RU/species_Mw[i_getCp];
        }
        
    }
    
    private:

    /**
     * Small helper class which provides wrapper to
     * Thermodynamics::speciesCpOverR() to be used with getTFromRhoE().
     */
    
    class Cp {
    public:
        Cp(const Thermodynamics& t, const double& T_vibration) : thermo(t) {
             Tv =  T_vibration;
//             p_cpt = new double [thermo.nSpecies()]; //DELETE IT
//             p_cpr = new double [thermo.nSpecies()]; //DELETE IT
//             p_cpv = new double [thermo.nSpecies()]; //DELETE IT
//             p_cpel = new double [thermo.nSpecies()]; //DELETE IT
        }
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(cp);
//          thermo.speciesCpOverR(T, Tv, T, Tv, Tv, cp, p_cpt, p_cpr, p_cpv, p_cpel);
//          thermo.speciesCpOverR(T, cp);
//          std::cout << "Cp_tot = " << cp[0]*RU/thermo.speciesMw(0) << "   " << cp[1]*RU/thermo.speciesMw(1) << endl;
//          std::cout << "Cp_tr = " << p_cpt[0]*RU/thermo.speciesMw(0) << "   " << p_cpt[1]*RU/thermo.speciesMw(1) << endl;
//          std::cout << "Cp_rot = " << p_cpr[0]*RU/thermo.speciesMw(0) << "   " << p_cpr[1]*RU/thermo.speciesMw(1) << endl;
//          std::cout << "Cp_vib = " << p_cpv[0]*RU/thermo.speciesMw(0) << "   " << p_cpv[1]*RU/thermo.speciesMw(1) << endl;
//          std::cout << "Cp_el = " << p_cpel[0]*RU/thermo.speciesMw(0) << "   " << p_cpel[1]*RU/thermo.speciesMw(1) << endl;
//          std::cout << endl;
        }

    private:
        const Thermodynamics& thermo;
        double Tv;
	double* p_cpt;
	double* p_cpr;
	double* p_cpv;
	double* p_cpel;
    }; // class Cp
    
    class Cpv {
    public:
        Cpv(const Thermodynamics& t) : thermo(t) {}
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(T, T, T, T, T, NULL, NULL, NULL, cp, NULL);
//	    std::cout << "For the vibrational loop = " << cp[0]*RU/thermo.speciesMw(0) << "    " << cp[1]*RU/thermo.speciesMw(0)<< endl;
        }
    private:
        const Thermodynamics& thermo;
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
	  //std::cout << "Total Enthalpy = " << h[0] << endl;
	  //T = 10000;
          thermo.speciesHOverRT(T, Tv, T, Tv, Tv, h, NULL, NULL, NULL, NULL, NULL);
//	  std::cout << "T = " << T << endl;
//	  std::cout << "Total Enthalpy = " << h[0] * 8.314/0.014 * 10000 << "    " << h[1] * 8.314/0.028 * 10000 << endl;
          //thermo.speciesHOverRT(T, h);
        }
    private:
        const Thermodynamics& thermo;
        double Tv;
    }; // class H
    
    class Hv {
    public:
        Hv(const Thermodynamics& t) : thermo(t) {
            p_enth = new double [t.nSpecies()];
	    p_trans = new double [t.nSpecies()];
        }
        ~Hv(){
            delete [] p_enth;
	    delete [] p_trans;
        }
        void operator () (double T, double* const h) const {
            thermo.speciesHOverRT(T, T, T, T, T, p_enth, p_trans, NULL, h, NULL, NULL);
        }
    private:
        const Thermodynamics& thermo;
        double* p_enth;
	double* p_trans;
    }; // class Hv

private:
    double* mp_work;
    Mutation::Transfer::TransferModel* mp_omega_VT;
    Mutation::Transfer::TransferModel* mp_omega_CV;
//    Mutation::Transfer::TransferModel* mp_omega_EV;
//    Mutation::Transfer::TransferModel* mp_omega_ET;
//    Mutation::Transfer::TransferModel* mp_omega_I;
    
}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqTTvStateModel, StateModel> chem_non_eq_TTv("ChemNonEqTTv");

    } // namespace Thermodynamics
} // namespace Mutation
