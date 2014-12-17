/**
 * @file StateModel.h
 *
 * @brief Definition of the StateModel class.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef THERMO_STATE_MODEL_H
#define THERMO_STATE_MODEL_H

#include "Constants.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

namespace Mutation {

class Mixture;

	namespace Transfer { class TransferModel; }

    namespace Thermodynamics {

/**
 * @defgroup statemodels State Models
 * List of available State Models.
 * @{
 */

/**
 * Base class for all state models.
 */
class StateModel
{
public:

    typedef Mutation::Mixture& ARGS;
    
    /**
     * Constructor initializes state.
     *
     * @param ns - number of species
     */
    StateModel(ARGS mix, const int nenergy, const int nmass);
    
    /**
     * Destructor.
     */
    virtual ~StateModel();
    
    /**
     * Returns the number of energy equations represented by this StateModel.
     */
    inline int nEnergyEqns() const { return m_nenergy; }

    /**
     * Returns the number of mass equations represented by this StateModel.
     */
    inline int nMassEqns() const { return m_nmass; }

    /**
     * Abstract method to set the state of the mixture which is dependent on the
     * type of the concrete StateModel class.
     *
     * @param p_mass - The "mass" vector (depends on concrete class)
     * @param p_energy - The "energy" vector (depends on concrete class)
     * @param vars - Index representing which variable set is given in the mass
     * and energy vectors
     */
    virtual void setState(
        const double* const p_mass, const double* const p_energy,
        const int vars = 0) = 0;
    
    /**
     * Returns the mixture translational temperature.
     */
    inline double T() const {
        return m_T;
    }
    
    /**
     * Returns the mixture vibrational temperature.
     */
    inline double Tv() const {
        return m_Tv;
    }
    
    /**
     * Returns the mixture electron temperature.
     */
    inline double Te() const {
        return m_Te;
    }
    
    /**
     * Returns the mixture rotational temperature.
     */
    inline double Tr() const {
        return m_Tr;
    }
    
    /**
     * Returns the mixture electronic temperature.
     */
    inline double Tel() const {
        return m_Tel;
    }
    
    /**
     * Returns the mixture static pressure.
     */
    inline double P() const {
        return m_P;
    }
    
    /**
     * Fills an array of the temperatures represented by this StateModel.
     */
    virtual void getTemperatures(double* const p_T) const {
        std::cerr << "getTemperatures()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * energy per unit mass.  The first n_species values correspond to the total energy
	 * vector.  The remaining n_species vectors correspond to the additional energy modes.
     */
    virtual void getEnergiesMass(double* const p_e) {
        std::cerr << "getEnergiesMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }
    
    virtual void getMixtureEnergiesMass(double* const p_e) {
        std::cerr << "getMixtureEnergiesMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * enthalpy per unit mass.  The first n_species values correspond to the total enthalpy
	 * vector.  The remaining n_species vectors correspond to the additional energy modes.
     */
    virtual void getEnthalpiesMass(double* const p_h) {
        std::cerr << "getEnthalpiesMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }
    
    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * cp per unit mass.  Each n_species vector corresponds to a temperature in the state
     * model. The first one is associated with the heavy particle translational temperature. 
     */
    
    virtual void getCpsMass(double* const p_Cp){
        std::cerr << "getCpsMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

    /**
     * Returns a vector of length n_species times n_energies with each corresponding
	 * cp per unit mass.  Each n_species vector corresponds to a temperature in the state
     * model. The first one is associated with the heavy particle translational temperature. 
     */
    
    virtual void getCvsMass(double* const p_Cp){
        std::cerr << "getCvsMass()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

	/**
	 * Assigns a unique temperature to each energy mode.
     */
    virtual void getTagModes(int* const p_tag) {
        std::cerr << "getTagModes()"
                  << " not implemented by this StateModel!" << std::endl;
        std::exit(1);
    }

    /**
     * Returns the species mole fractions.
     */
    inline const double* const X() const {
        return mp_X;
    }
    
    /**
     * This function provides the total energy transfer source terms
     */
    virtual void energyTransferSource(double* const p_omega);
    
protected:

    void addTransferTerm(int i, const std::string& model);

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
     * @param p_work    work array, at least number of species long
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
        double* const p_work,
        const double alpha = 0.0,
        const double atol = 1.0e-12,
        const double rtol = 1.0e-12,
        const int max_iters = 100)
    {
        const double rhoe_over_Ru = rhoe/RU;
        const double tol = rtol*std::abs(rhoe_over_Ru) + atol;

        double f, fp, dT;

        // Compute initial value of f
        h(T, p_work);
        f = alpha;
        for (int i = 0; i < m_ns; ++i)
            f += mp_X[i]*p_work[i];
        f = T*f - rhoe_over_Ru;

        int iter = 0;
        //cout << iter << " " << f << " " << T << endl;
        while (std::abs(f) > tol) {
            // Check for max iterations
            if (iter++ == max_iters) {
                using std::cerr;
                std::cerr << "Exceeded max iterations when computing temperature!\n";
                std::cerr << "res = " << f / rhoe_over_Ru << ", T = " << T << std::endl;
                return false;
            }

            // Compute df/dT
            cp(T, p_work);
            fp = alpha;
            for (int i = 0; i < m_ns; ++i)
                fp += mp_X[i]*p_work[i];

            // Update T
            dT = f/fp;
            while (dT >= T) dT *= 0.5; // prevent non-positive T
            T -= dT;

            // Recompute f
            h(T, p_work);
            f = alpha;
            for (int i = 0; i < m_ns; ++i)
                f += mp_X[i]*p_work[i];
            f = T*f - rhoe_over_Ru;
            //cout << iter << " " << f << " " << T << endl;
        }

        // Let the user know if we converged or not
        return true;
    }

protected:

    Mutation::Mixture& m_mix;
    const int m_ns;
    const int m_nenergy;
    const int m_nmass;
    
    double m_T;
    double m_Tv;
    double m_Tr;
    double m_Tel;
    double m_Te;
    double m_P;
    
    double* mp_X;
    

    std::vector< std::pair<int, Mutation::Transfer::TransferModel*> >
        m_transfer_models;

}; // class StateModel

/// @}

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_STATE_MODEL_H

