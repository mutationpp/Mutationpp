/**
 * @file ChemNonEqStateModel.h
 *
 * @brief Provides ChemNonEqStateModel.
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
     *   2: primitive set 2 (species mass fractions, {P, T} array)
     */
    void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0)
    {
        const int ns = m_thermo.nSpecies();

        // Compute the species concentrations which are used throughout this
        // method regardless of variable set.  Note that in the case where
        // the p_mass vector is mass fractions, this is just dividing by the
        // molecular weights as the first step in converting to mole fractions.
        double conc = 0.0;
        for (int i = 0; i < ns; ++i) {
            // Check that species densities are at least positive
            assert(p_mass[i] >= 0.0);
            mp_X[i] = std::max(p_mass[i] / m_thermo.speciesMw(i), 0.0);
            conc += mp_X[i];
        }

        // Compute the temperature and pressure and make sure the variable set
        // is implemented
        switch (vars) {
        case 0:
            // Solve energy equation for temperature
            getTFromRhoE(
                Cp(m_thermo), H(m_thermo), p_energy[0], m_T, mp_work, -conc);
            m_P = RU * m_T * conc;
            break;

        case 1:
            // Check that temperature is at least positive
            assert(p_energy[0] > 0.0 && "Temperature is negative!");
            m_T = p_energy[0];
            m_P = RU * m_T * conc;
            break;

        case 2:
            // Check temperature and pressure are positive
            assert(p_energy[0] > 0.0 && "Pressure is negative!");
            assert(p_energy[1] > 0.0 && "Temperature is negative!");
            m_P = p_energy[0];
            m_T = p_energy[1];
            break;

        default:
            throw InvalidInputError("variable set", vars)
                << "This variable-set is not implemented in ChemNonEqStateModel"
                << ". Possible variable-sets are:\n"
                << "  0: (species densities, static energy density)\n"
                << "  1: (species densities, temperature)\n"
                << "  2: (species mass fractions, pressure and temperature)";
        }

        // All other temperatures are the same
        m_Tr = m_Tv = m_Tel = m_Te = m_T;

        // Compute the species mole fractions
        for (int i = 0; i < ns; ++i)
            mp_X[i] /= conc;
    }

    void getTemperatures(double* const p_T) const {
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
