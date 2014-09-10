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
        for (int i = 0; i < ns; ++i) {
            // Check that species densities are at least positive
            assert(p_mass[i] >= 0.0);
            mp_X[i] = p_mass[i] / m_thermo.speciesMw(i);
        }

        // Compute the temperature and make sure the variable set is implemented
        switch (vars) {
        case 0: {
            // Solve nonlinear system to get T from rho*e using Newton's method
            // Use last state update as initial guess for T
            const int max_iters = 50;
            const double tol = 1.0e-12;
            double f, fp, dT;

            m_thermo.speciesHOverRT(m_T, mp_work);
            f = 0.0;
            for (int i = 0; i < ns; ++i)
                f += mp_X[i]*(mp_work[i] - 1.0);
            f = RU*m_T*f - *p_energy;

            int iter = 0;
            while (std::abs(f) > tol * m_T) {
                // Print warning if this is taking too long
                if (++iter % max_iters == 0) {
                    cout << "setState() taking too many iterations for ChemNonEq1T StateModel!"
                         << " It is likely that the input arguments are not feasible..." << endl;
                    cout << "Species densities [kg/m^3]:" << endl;
                    for (int i = 0; i < ns; ++i)
                        cout << "  " << setw(20) << m_thermo.speciesName(i) << " " << p_mass[i] << endl;
                    cout << "Energy density [J/m^3]:" << endl;
                    cout << p_energy[0] << endl;
                }

                // Compute df/dT
                m_thermo.speciesCpOverR(m_T, mp_work);
                fp = 0.0;
                for (int i = 0; i < ns; ++i)
                    fp += mp_X[i]*(mp_work[i] - 1.0);
                fp *= RU;

                // Update T
                dT = f/fp;
                while (dT > m_T) dT *= 0.5; // prevent negative T
                m_T -= dT;

                // Recompute f
                m_thermo.speciesHOverRT(m_T, mp_work);
                f = 0.0;
                for (int i = 0; i < ns; ++i)
                    f += mp_X[i]*(mp_work[i] - 1.0);
                f = RU*m_T*f - *p_energy;
            }

            break;
        }
        case 1:
            // Check that temperature is at least positive
            assert(*p_energy > 0.0);
            m_T = *p_energy;
            break;
        default:
            cout << "Variable-set " << vars << " not implemented in StateModel!" << endl;
            exit(1);
        }

        // All other temperatures are the same
        m_Tr = m_Tv = m_Tel = m_Te = m_T;

        // Compute the pressure and species mole fractions from T and rho_i
        m_P = 0.0;
        for (int i = 0; i < ns; ++i)
            m_P += mp_X[i];

        for (int i = 0; i < ns; ++i)
            mp_X[i] /= m_P;
        m_P *= RU * m_T;
    }

private:

    double* mp_work;

}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqStateModel, StateModel> chem_non_eq_1T("ChemNonEq1T");

    } // namespace Thermodynamics
} // namespace Mutation
