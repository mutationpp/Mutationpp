#include "StateModel.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * Chemical Non-equilibrium 2-T StateModel.  Thermal nonequilibrium with
 * T = Tr and Tv = Tel = Te with finite rate chemistry.
 */
class ChemNonEqTTvStateModel : public StateModel
{
public:

    ChemNonEqTTvStateModel(const Thermodynamics& thermo)
        : StateModel(thermo, 2, thermo.nSpecies())
    {
        // Allocate work storage in a single block
        mp_work1 = new double [2*thermo.nSpecies()];
        mp_work2 = mp_work1 + thermo.nSpecies();
    }

    ~ChemNonEqTTvStateModel()
    {
        delete [] mp_work1;
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
        const int offset = ns - m_thermo.nHeavy();
        const double tol = 1.0e-8;

        // Compute the species concentrations which are used through out this
        // method regardless of variable set (store in mole fraction array
        // because it is efficient at the end to avoid a copy and does not
        // require another work array).
        for (int i = 0; i < ns; ++i)
            mp_X[i] = p_mass[i] / m_thermo.speciesMw(i);

        // Compute the temperatures and make sure the requested variable set is
        // actually implemented
        switch (vars) {
        case 0: {
            // First step is to solve one nonlinear equation to get Tv from the
            // vibrational energy density equation
            double fv, fvp;
            ThermoDB* p_thermodb = m_thermo.thermoDB();

            p_thermodb->hv(m_Tv, mp_work1);
            fv = 0.0;
            for (int i = offset; i < ns; ++i)
                fv += mp_X[i]*mp_work1[i];
            fv = RU*m_Tv*fv - p_energy[1];

            while (std::abs(fv) > tol) {
                // Compute df/dT
                p_thermodb->cpv(m_Tv, mp_work1);
                fvp = 0.0;
                for (int i = offset; i < ns; ++i)
                    fvp += mp_X[i]*mp_work1[i];
                fvp *= RU;

                // Update T
                m_Tv -= fv/fvp;
                if (m_Tv <= 0.0) m_Tv = 1.0;

                // Recompute f
                p_thermodb->hv(m_Tv, mp_work1);
                fv = 0.0;
                for (int i = offset; i < ns; ++i)
                    fv += mp_X[i]*mp_work1[i];
                fv = RU*m_Tv*fv - p_energy[1];
            }

            // Next compute the (trans.-rot.) energy density
            double ert_density = 0.0;
            p_thermodb->hel(m_Tv, mp_work2);
            for (int i = offset; i < ns; ++i)
                ert_density += mp_X[i]*(mp_work1[i]+mp_work2[i]);
            if (offset == 1)
                ert_density += 1.5*mp_X[0];
            ert_density = p_energy[0] - RU*m_Tv*ert_density;




            m_T = m_Tv;

            break;
        }
        case 1:
            m_T  = p_energy[0];
            m_Tv = p_energy[1];
            break;
        default:
            cout << "Variable set not implemented in StateModel!" << endl;
            exit(1);
        }

        // Set Tr, Tel, and Te
        m_Tr = m_T;
        m_Tel = m_Te = m_Tv;

        // Compute the pressure and species mole fractions from T and rho_i
        m_P = 0.0;
        for (int i = 0; i < ns; ++i)
            m_P += mp_X[i];

        for (int i = 0; i < ns; ++i)
            mp_X[i] /= m_P;
        m_P *= RU * m_T;
    }

private:

    double* mp_work1;
    double* mp_work2;

}; // class ChemNonEqStateModel

// Register the state model
Utilities::Config::ObjectProvider<
    ChemNonEqTTvStateModel, StateModel> chem_non_eq_TTv("ChemNonEqTTv");

    } // namespace Thermodynamics
} // namespace Mutation
