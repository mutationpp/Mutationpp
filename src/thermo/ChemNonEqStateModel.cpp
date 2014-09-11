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

    // Define the H and Cp providers
    class Cp {
    public:
        const Thermodynamics& thermo;
        Cp(const Thermodynamics& t) : thermo(t) {}
        void operator () (double T, double* const cp) const {
            thermo.speciesCpOverR(T, cp);
        }
    };

    class H {
    public:
        const Thermodynamics& thermo;
        H(const Thermodynamics& t) : thermo(t) {}
        void operator () (double T, double* const cp) const {
            thermo.speciesHOverRT(T, cp);
        }
    };

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
            getTFromRhoE(Cp(m_thermo), H(m_thermo), p_energy[0], m_T, -conc, 0.0);
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

protected:

    /**
     * Solves the general form of an energy equation
     * \f[ f(T) = T \left[\sum_i \tilde{\rho}_i \left(\frac{H_i(T)}{R_uT}\right)
     *            + \alpha\right] - \frac{\rho e}{R_u} = 0 \f]
     * for the temperature given the \f$C_{p,i}/R_u\f$ and \f$H_i/R_uT\f$
     * function providers using a Newton-Rhapson iterative procedure.  This
     * function assumes that the mp_X array is filled with the species molar
     * densities (concentrations) prior to being called. The convergence
     * criteria for the Newton iterations is taken to be
     * \f[ |f(T)| < \epsilon_r\left|\frac{\rho e}{R_u}\right| + \epsilon_a \f]
     * where \f$\epsilon_a\f$ and \f$\epsilon_r\f$ are absolute and relative
     * tolerances respectively.  Note that if the total energy equation is
     * being solved, \f$\alpha\f$ should be set to the negative mixture
     * concentration, \f$-\sum \tilde{\rho}_i\f$. For internal energy equations,
     * \f$\alpha = 0\f$.
     *
     * @param cp        class providing species specific heats
     * @param h         class providing species enthalpies
     * @param rhoe      value of \f$\rho e\f$
     * @param T         initial guess on input, output is the solution
     * @param alpha     arbitrary temperature coefficient in energy equation
     * @param atol      absolute tolerance, \f$\epsilon_a\f$, on \f$f\f$
     * @param rtol      relative tolerance, \f$\epsilon_r\f$, on \f$f\f$
     * @param max_iters maximum number of iterations (warns user if exceeded)
     *
     * @return false if max iterations exceeded, true otherwise
     */
    template <typename CpProvider, typename HProvider>
    bool getTFromRhoE(
        const CpProvider& cp,
        const HProvider& h,
        const double rhoe,
        double& T,
        const double alpha = 0.0,
        const double atol = 1.0e-12,
        const double rtol = 1.0e-12,
        const int max_iters = 10)
    {
        const int ns = m_thermo.nSpecies();
        const double rhoe_over_Ru = rhoe/RU;
        const double tol = rtol*std::abs(rhoe_over_Ru) + atol;

        double f, fp, dT;

        // Compute initial value of f
        h(T, mp_work);
        f = alpha;
        for (int i = 0; i < ns; ++i)
            f += mp_X[i]*mp_work[i];
        f = T*f - rhoe_over_Ru;

        int iter = 0;
        cout << iter << " " << f << " " << T << endl;
        while (std::abs(f) > tol) {
            // Check for max iterations
            if (iter++ == max_iters) {
                using std::cerr;
                cerr << "Exceeded max iterations when computing temperature!\n";
                cerr << "res = " << f / rhoe_over_Ru << ", T = " << T << endl;
                return false;
            }

            // Compute df/dT
            cp(T, mp_work);
            fp = alpha;
            for (int i = 0; i < ns; ++i)
                fp += mp_X[i]*mp_work[i];

            // Update T
            dT = f/fp;
            while (dT >= T) dT *= 0.5; // prevent non-positive T
            T -= dT;

            // Recompute f
            h(T, mp_work);
            f = alpha;
            for (int i = 0; i < ns; ++i)
                f += mp_X[i]*mp_work[i];
            f = T*f - rhoe_over_Ru;

            cout << iter << " " << f << " " << T << endl;
        }

        // Let the user know if we converged or not
        return true;
    }

private:

    double* mp_work;

}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqStateModel, StateModel> chem_non_eq_1T("ChemNonEq1T");

    } // namespace Thermodynamics
} // namespace Mutation
