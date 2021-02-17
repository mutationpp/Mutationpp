/**
 * @file EquilStateModel.h
 *
 * @brief Provides EquilStateModel class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "StateModel.h"
#include "Thermodynamics.h"
#include "Transport.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * @ingroup statemodels
 *
 * StateModel for thermochemical equilibrium.
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
        mp_work = new double [ns];
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
        delete [] mp_work;
    }

    /**
     * Sets the mixture state by computing the equilibrium composition at the
     * given mixture conditions based on the given variable set. Variable sets
     * can be the following:
     *   0: conserved variables (mixture density, static energy density)
     *   1: primitive set 1 (pressure, temperature)
     *   2: primitive set 2 (elemental mole fractions,  {P, T} array)
     */
    virtual void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0)
    {
        // Determine temperature, pressure, and mole fractions depending on
        // variable set
        switch (vars) {
            // Given density and energy density (conserved variables)
            case 0: {
                assert(p_mass[0]   > 0.0);
                assert(p_energy[0] > 0.0);

                const int ns = m_thermo.nSpecies();
                const int max_iters = 50;
                const double tol = 1.0e-12;
                const double tolAbs = 1.0e-10;
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
                while ((std::max(std::abs(f1/(std::abs(rhoe)+tolAbs)),std::abs(f2/(rho+tolAbs))) > tol)
                        && (std::max(std::abs(f1),std::abs(f2)) > tolAbs)) {
                    // Print warning if this is taking too long
                    if (++iter % max_iters == 0) {
                        std::cout << "setState() taking too many iterations for Equil StateModel!"
                             << " It is likely that the input arguments are not feasible..." << std::endl;
                        std::cout << "density [kg/m^3] = " << rho << ", energy [J/m^3] = " << rhoe << std::endl;
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
                    dT = -f1/dfdt;
                    while (dT > m_T) dT *= 0.5; // prevent negative T
                    m_T -= dT;

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

            // Given pressure and temperature (using default elemental fractions)
            case 1:
                assert(p_mass[0]   > 0.0);
                assert(p_energy[0] > 0.0);

                m_T = p_energy[0];
                m_P = p_mass[0];
                m_thermo.equilibriumComposition(m_T, m_P, mp_X);
                break;

            // Given elemental mole fractions and temperature and pressure
            case 2:
                assert(p_energy[0] > 0.0);
                assert(p_energy[1] > 0.0);

                m_T = p_energy[1];
                m_P = p_energy[0];
                m_thermo.equilibriumComposition(m_T, m_P, p_mass, mp_X);
                break;

            // Unknown variable set
            default:
                throw InvalidInputError("variable set", vars)
                    << "This variable-set is not implemented in EquilStateModel"
                    << ". Possible variable-sets are:\n"
                    << "  0: (mixture density, static energy density)\n"
                    << "  1: (pressure, temperature)\n"
                    << "  2: (element mole fractions, pressure and temperature)";
        }

        // Set the remaining temperatures
        m_Tr = m_Tv = m_Tel = m_Te = m_T;
    }

    void getTemperatures(double* const p_T) const
    {
        p_T[0] = m_T;
    }

    void getEnergiesMass(double* const p_e)
	{
		const int ns = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work);

        for(int i = 0; i < ns; ++i)
            p_e[i] = (mp_work[i]  - 1.0)*m_T*RU/m_thermo.speciesMw(i);
    }

    void getEnthalpiesMass(double* const p_h)
    {
		const int ns = m_thermo.nSpecies();
        m_thermo.speciesHOverRT(mp_work);

        for(int i = 0; i < ns; ++i)
            p_h[i] = mp_work[i]*m_T*RU/m_thermo.speciesMw(i);
    }

    void getCpsMass(double* const p_Cp)
    {
        const int ns = m_thermo.nSpecies();
        m_thermo.speciesCpOverR(m_T, mp_work);

        for(int i = 0; i < ns; ++i)
            p_Cp[i] = mp_work[i]*RU/m_thermo.speciesMw(i);
    }

	void getCvsMass(double* const p_Cv)
	{
        const int ns = m_thermo.nSpecies();
        m_thermo.speciesCpOverR(m_T, mp_work);

        for(int i = 0; i < ns; ++i)
            p_Cv[i] = (mp_work[i]-1.0)*RU/m_thermo.speciesMw(i);
    }


private:

    double* mp_h;
    double* mp_work;
    double* mp_cp;
    double* mp_dxdt;
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilStateModel, StateModel> equil_sm("Equil");

/**
 * @ingroup statemodels
 *
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
          //assert(vars == 1);
          EquilStateModel::setState(p_P, p_T, 1);
     }
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilTPStateModel, StateModel> equiltp_sm("EquilTP");


    } // namespace Thermodynamics
} // namespace Mutation
