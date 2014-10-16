#include "StateModel.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * Chemical Non-equilibrium StateModel.  Mixture is in thermal equilibrium but
 * chemistry source terms are finite-rate.
 */
class ChemNonEqStateModel : public StateModel
{
public:

    ChemNonEqStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 1, thermo.nSpecies())
    {
        mp_work = new double [thermo.nSpecies()];
    }

    ~ChemNonEqStateModel()
    {
        delete [] mp_work;
    }

    /**
     * Sets the current state of the mixture.  Variable sets can be the
     * following:
     *   0: conserved variables (species densities, total energy density)
     *   1: primitive set 1 (species densities, mixture temperature)
     */
    void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0)
    {
        const int ns = m_thermo.nSpecies();

        // Compute the species concentrations which are used throughout this
        // method regardless of variable set
        double conc = 0.0;
        for (int i = 0; i < ns; ++i) {
            // Check that species densities are at least positive
            assert(p_mass[i] >= 0.0);
            mp_X[i] = p_mass[i] / m_thermo.speciesMw(i);
            conc += mp_X[i];
        }

        // Compute the temperature and make sure the variable set is implemented
        switch (vars) {
        case 0:
            // Solve energy equation
            getTFromRhoE(
                Cp(m_thermo), H(m_thermo), p_energy[0], m_T, mp_work, -conc);
            break;
        case 1:
            // Check that temperature is at least positive
            assert(p_energy[0] > 0.0);
            m_T = p_energy[0];
            break;
        default:
            cout << "Variable-set " << vars << " not implemented in StateModel!" << endl;
            exit(1);
        }

        // All other temperatures are the same
        m_Tr = m_Tv = m_Tel = m_Te = m_T;

        // Compute the pressure and species mole fractions from T and rho_i
        for (int i = 0; i < ns; ++i)
            mp_X[i] /= conc;
        m_P = RU * m_T * conc;
    }

    void getEnergyMass(double* const p_e) {

        int n_species = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work);
        double species_Mw[n_species];

        for(int i_energy_mass = 0; i_energy_mass < m_thermo.nSpecies(); ++i_energy_mass){
            species_Mw[i_energy_mass] = m_thermo.speciesMw(i_energy_mass);
            p_e[i_energy_mass] = mp_work[i_energy_mass] * m_T * RU /species_Mw[i_energy_mass] - RU/species_Mw[i_energy_mass] *m_T;
        }

    }

    void getEnthalpyMass(double* const p_h) {

        int n_species = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work);
        double species_Mw[n_species];

        for(int i_enthalpy_mass = 0; i_enthalpy_mass < n_species; ++i_enthalpy_mass){
            species_Mw[i_enthalpy_mass] = m_thermo.speciesMw(i_enthalpy_mass);
            p_h[i_enthalpy_mass] = mp_work[i_enthalpy_mass] * m_T * RU /species_Mw[i_enthalpy_mass];
        }

    }

    void getCpMass(double* const p_Cp) {

        int n_species = m_thermo.nSpecies();
        m_thermo.speciesCpOverR(m_T, mp_work);
        double species_Mw[n_species];

        for(int i_getCp = 0; i_getCp < n_species; ++i_getCp){
            species_Mw[i_getCp] = m_thermo.speciesMw(i_getCp);
            p_Cp[i_getCp] = mp_work[i_getCp] *RU/species_Mw[i_getCp];
        }

    }

private:

    /**
     * Small helper class which provides wrapper to
     * Thermodynamics::speciesCpOverR() to be used with getTFromRhoE().
     */
    class Cp {
    public:
        Cp(const Thermodynamics& t) : thermo(t) {}
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(T, cp);
        }
    private:
        const Thermodynamics& thermo;
    };


    /**
     * Small helper class which provides wrapper to
     * Thermodynamics::speciesHOverRT() to be used with getTFromRhoE().
     */
    class H {
    public:
        H(const Thermodynamics& t) : thermo(t) {}
        void operator () (double T, double* const h) const {
            thermo.speciesHOverRT(T, h);
        }
    private:
        const Thermodynamics& thermo;
    };

private:

    double* mp_work;

}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqStateModel, StateModel> chem_non_eq_1T("ChemNonEq1T");

    } // namespace Thermodynamics
} // namespace Mutation
