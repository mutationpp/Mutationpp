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
  
     Mutation::Transfer::TransferModel* p_transfer_model;

    ChemNonEqTTvStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 2, thermo.nSpecies())
    {
        mp_work = new double [thermo.nSpecies()];
    }

    ~ChemNonEqTTvStateModel()
    {
        delete [] mp_work;
        delete p_transfer_model;
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
                Cpv(m_thermo), Hv(m_thermo), p_energy[1], m_Tv, mp_work, 0.0);

            // Next compute the temperature using the total energy density

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

        p_transfer_model = new Mutation::Transfer::TransferModelVT(m_thermo); // I include a copy of thermo
        double output;
        p_transfer_model->source(&output);
	
	std::cout << "Source VT = " << output << endl;

    }

      Mutation::Transfer::TransferModel* createTransferModel(
      Thermodynamics& thermo,
      Mutation::Transport::Transport& transport,
      Mutation::Kinetics::Kinetics& kinetics)
    {
      p_transfer_model = new Mutation::Transfer::TransferModelVT(thermo);  // DO NOT FORGET TO DELETE IT
      return p_transfer_model;
      
    }
    
    void getTemperatures(double* const p_T) const {
        p_T[0] = T();
        p_T[1] = Tv();
    }
    
    void getEnergyDensities(double* const p_rhoe) {
        std::cerr << "getEnergyDensities()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
/*        double rho = m_thermo.density();  
        int n_species = m_thermo.nSpecies();

        m_thermo.speciesHOverRT(p_rhoe, NULL, NULL, p_rhoe+n_species, NULL, NULL);

        const double* const Y = m_thermo.Y();

        for(int i_energy_densities = 0; i_energy_densities < m_thermo.nSpecies(); ++i_energy_densities){
            p_rhoe[i_energy_densities] *=  m_T * RU * rho * Y[i_energy_densities];
            p_rhoe[i_energy_densities+n_species] *=  m_T * RU * rho * Y[i_energy_densities];
        }
*/
    }
    
    void getEnthalpyDensities(double* const p_rhoh) {
        
        double rho = m_thermo.density();  
        int n_species = m_thermo.nSpecies();

        m_thermo.speciesHOverRT(p_rhoh, NULL, NULL, p_rhoh+n_species, NULL, NULL);

        const double* const Y = m_thermo.Y();

        for(int i_enthalpy_densities = 0; i_enthalpy_densities < m_thermo.nSpecies(); ++i_enthalpy_densities){
            
            p_rhoh[i_enthalpy_densities] *=  m_T * RU * rho * Y[i_enthalpy_densities];
            p_rhoh[i_enthalpy_densities+n_species] *=  m_T * RU * rho * Y[i_enthalpy_densities];
        }

    }
    
    void getCp(double* const p_Cp) {
        int n_species = m_thermo.nSpecies();
        m_thermo.speciesCpOverR(m_T, m_Te, m_Tr, m_Tv, m_Tel, p_Cp, NULL, NULL, p_Cp+n_species, NULL);

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
        }
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(T, Tv, T, Tv, Tv, cp, NULL, NULL, NULL, NULL);
        }

    private:
        const Thermodynamics& thermo;
        double Tv;
    }; // class Cp
    
    class Cpv {
    public:
        Cpv(const Thermodynamics& t) : thermo(t) {}
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(T, T, T, T, T, NULL, NULL, NULL, cp, NULL);
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
            thermo.speciesHOverRT(T, Tv, T, Tv, Tv, h, NULL, NULL, NULL, NULL, NULL);
        }
    private:
        const Thermodynamics& thermo;
        double Tv;
    }; // class H
    
    class Hv {
    public:
        Hv(const Thermodynamics& t) : thermo(t) {
            p_enth = new double [t.nSpecies()];
        }
        ~Hv(){
            delete [] p_enth;
        }
        void operator () (double T, double* const h) const {
            thermo.speciesHOverRT(T, T, T, T, T, p_enth, NULL, NULL, h, NULL, NULL);
        }
    private:
        const Thermodynamics& thermo;
        double* p_enth;
    }; // class Hv

private:
    double* mp_work;

}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqTTvStateModel, StateModel> chem_non_eq_TTv("ChemNonEqTTv");

    } // namespace Thermodynamics
} // namespace Mutation
