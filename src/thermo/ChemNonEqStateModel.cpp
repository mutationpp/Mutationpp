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
        const double atol = 1.0e-12;
        const double rtol = 1.0e-12;
        const int max_iters = 50;

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
        case 0: {
            // Solve nonlinear system to get T from rho*e using Newton's method
            // Use last state update as initial guess for T
            const double rhoe_over_Ru = p_energy[0]/RU;
            double f, fp, dT;
            
            m_thermo.speciesHOverRT(m_T, mp_work);
            f = -conc;
            for (int i = 0; i < ns; ++i)
                f += mp_X[i]*mp_work[i];
            f = m_T*f - rhoe_over_Ru;

            int iter = 0;
            const double tol = rtol*std::abs(rhoe_over_Ru) + atol;
            while (std::abs(f) > tol) {
                // Print warning if this is taking too long
                if (++iter == max_iters) {
                    cout << "Reached maximum iterations in ChemNonEq1T::setState()!" << endl;
                    cout << "f = " << f / rhoe_over_Ru << ", T = " << m_T << endl;
                    break;
                }

                // Compute df/dT
                m_thermo.speciesCpOverR(m_T, mp_work);
                fp = -conc;
                for (int i = 0; i < ns; ++i)
                    fp += mp_X[i]*mp_work[i];

                // Update T
                dT = f/fp;
                while (dT >= m_T) dT *= 0.5; // prevent nonpositive T
                m_T -= dT;

                // Recompute f
                m_thermo.speciesHOverRT(m_T, mp_work);
                f = -conc;
                for (int i = 0; i < ns; ++i)
                    f += mp_X[i]*mp_work[i];
                f = m_T*f - rhoe_over_Ru;
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
        for (int i = 0; i < ns; ++i)
            mp_X[i] /= conc;
        m_P = RU * m_T * conc;
    }

private:

    double* mp_work;

}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqStateModel, StateModel> chem_non_eq_1T("ChemNonEq1T");

    } // namespace Thermodynamics
} // namespace Mutation
