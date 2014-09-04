#include "StateModel.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * StateModel for thermochemical equilibrium.
 * @author J.B. Scoggins
 * @author Pierre Schrooyen
 */
class EquilStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    EquilStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 1, 1)
    {
        const int ns = m_thermo.nSpecies();

        // Allocate work storage
        mp_h = new double [ns*3];
        mp_cp = mp_h + ns;
        mp_dxdt = mp_cp + ns;

        // Initial solution
        double T = 300.0;
        double P = ONEATM;
        setState(&P, &T, 1);
    }

    virtual ~EquilStateModel()
    {
        delete [] mp_h; // other arrays allocated by this one
    }

    /**
     * Sets the mixture state by computing the equilibrium composition at the
     * given mixture conditions based on the given variable set. Variable sets
     * can be the following:
     *   0: conserved variables (mixture density, static energy density)
     *   1: primitive set 1 (pressure, temperature)
     */
    virtual void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0)
    {
        assert(p_mass[0] > 0.0);
        assert(p_energy[0] > 0.0);

        // Determine temperature, pressure, and mole fractions depending on
        // variable set
        switch (vars) {
            // Given density and energy density (conserved variables)
            case 0: {
                const int ns = m_thermo.nSpecies();
                const int max_iters = 50;
                const double tol = 1.0e-12;
                const double rho  = p_mass[0];
                const double rhoe = p_energy[0];

                double mw, h, cpeq, dfdt, dmwdt, f1, f2, dT;

                // Use the previous solution as the initial guess
                //m_thermo.equilibriumComposition(m_T, m_P, mp_X);
                m_thermo.speciesHOverRT(m_T, mp_h);
                mw = m_thermo.mixtureMw();

                h = 0.0;
                for (int i = 0; i < ns; ++i)
                    h += mp_h[i]*mp_X[i];
                h *= RU*m_T;
                f1 = rhoe - rho*h/mw + m_P;
                f2 = rho - m_P*mw/(RU*m_T);

                int iter = 0;
                while (std::max(std::abs(f1/rhoe),std::abs(f2/rho)) > tol) {
                    // Print warning if this is taking too long
                    if (iter++ % max_iters == 0) {
                        cout << "setState() taking too many iterations for Equil StateModel!"
                             << " It is likely that the input arguments are not feasible..." << endl;
                        cout << "density [kg/m^3] = " << rho << ", energy [J/m^3] = " << rhoe << endl;
                    }

                    // Compute df1/dT which is the main term
                    m_thermo.speciesCpOverR(m_T, mp_cp);

                    for (int i = 0; i < ns; ++i)
                        mp_dxdt[i] = -mp_h[i] / m_T;
                    m_thermo.equilSolver()->dXdg(mp_dxdt, mp_dxdt);

                    cpeq = 0;
                    for (int i = 0; i < ns; ++i)
                        cpeq += mp_h[i]*mp_dxdt[i];
                    cpeq *= m_T;
                    for (int i = 0; i < ns; ++i)
                        cpeq += mp_cp[i]*mp_X[i];
                    cpeq *= RU;

                    dmwdt = 0.0;
                    for (int i = 0; i < ns; ++i)
                        dmwdt += m_thermo.speciesMw(i)*mp_dxdt[i];

                    dfdt = rho*(cpeq*mw - h*dmwdt)/(mw*mw);

                    // Update T
                    dT = f1/dfdt;
                    while (dT > m_T) dT *= 0.5; // prevent negative T
                    m_T = m_T + dT;

                    // Update P (lagging Mwmix)
                    m_P = rho*RU*m_T/mw;

                    // Update the composition
                    m_thermo.equilibriumComposition(m_T, m_P, mp_X);

                    // Recompute f (and mw, h)
                    m_thermo.speciesHOverRT(m_T, mp_h);
                    mw = m_thermo.mixtureMw();

                    h = 0.0;
                    for (int i = 0; i < ns; ++i)
                        h += mp_h[i]*mp_X[i];
                    h *= RU*m_T;
                    f1 = rhoe - rho*h/mw + m_P;
                    f2 = rho - m_P*mw/(RU*m_T);
                }

            }
            break;

            // Given pressure and temperature
            case 1:
                m_T = p_energy[0];
                m_P = p_mass[0];
                m_thermo.equilibriumComposition(m_T, m_P, mp_X);
                break;

            // Unknown variable set
            default:
                cout << "Unknown variable set in Equil StateModel!" << endl;
                exit(1);
        }

        // Set the remaining temperatures
        m_Tr = m_Tv = m_Tel = m_Te = m_T;
    }

private:

    double* mp_h;
    double* mp_cp;
    double* mp_dxdt;
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilStateModel, StateModel> equil_sm("Equil");


/**
 * This class is just provided for backward compatibility purposes and may be
 * removed in the future.  Use the EquilStateModel instead which is more general
 * and will continue to be maintained in the future.
 *
 * @see EquilStateModel
 */
class EquilTPStateModel : public EquilStateModel
{
public:
    EquilTPStateModel(const Thermodynamics& thermo)
        : EquilStateModel(thermo)
    { }

    /**
     * Sets the state of the equilibrium mixture given temperature and pressure.
     */
    virtual void setState(
        const double* const p_T, const double* const p_P, const int vars = 1)
    {
        assert(vars == 1);
        EquilStateModel::setState(p_P, p_T, 1);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilTPStateModel, StateModel> equiltp_sm("EquilTP");

    } // namespace Thermodynamics
} // namespace Mutation
