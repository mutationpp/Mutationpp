/**
 * @file Kinetics.h
 *
 * @brief Declaration of Kinetics class.
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

#ifndef KINETICS_H
#define KINETICS_H

#include "StoichiometryManager.h"
#include "ThirdBodyManager.h"
#include "RateManager.h"
#include "JacobianManager.h"
#include "Reaction.h"
#include "Thermodynamics.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Manages the computation of chemical source terms for an entire reaction
 * mechanism.
 *
 * @todo Currently, intermediate values are always recomputed even if they don't
 * need to be.  An detailed look should be made into how to save computations
 * when, for instance, the temperature does not change.  This also holds true
 * for the jacobian.  For example, if the source terms were just computed, there
 * is no need to recompute the reaction rate coefficients in the jacobian
 * evaluation.
 */
class Kinetics
{
public:

    /**
     * Constructor which takes a reference to a Thermodynamics object and the
     * the full file name path to the mechanism data file.
     */
    Kinetics(
        const Mutation::Thermodynamics::Thermodynamics& thermo, 
        std::string mechanism);
    
    /**
     * Destructor.
     */
    ~Kinetics();
    
    /**
     * Returns the number of reactions in the mechanism.
     */
    size_t nReactions() const {
        return m_reactions.size();
    }
    
    /**
     * Returns the vector of Reaction objects associated with this kinetics
     * manager.
     */
    const std::vector<Reaction>& reactions() const {
        return m_reactions;
    }
    
    /**
     * Computes the change in any given species quantity across each reaction.
     */
    void getReactionDelta(const double* const p_s, double* const p_r) const;
    
    /**
     * Fills the vector keq with the equilibrium constant for each reaction 
     * \f$ K_{C,j} = \left(\frac{P_{atm}}{R_u T}\right)^{\Delta \nu_j}
     * \exp\left(-\frac{\Delta {G^\circ}_j}{R_u T}\right) \f$.
     *
     * @param T    the temperature in K
     * @param keq  on return, \f$ K_{C,j}(T) \f$
     */
    //void equilibriumConstants(const double T, double* const p_keq);
    
    /**
     * Fills the vector kf with the forward rate coefficients \f$ k_{f,j} \f$
     * for each reaction.  The coefficients are computed based on the rate laws
     * assigned to each individual reaction.
     *
     * @param kf  on return, \f$ k_{f,j} \f$
     */
    void forwardRateCoefficients(double* const p_kf);
    
    /**
     * Fills the vector kb with the backward rate coefficients 
     * \f$ k_{b,j} = k_{f,j} / K_{C,j} \f$ for each reaction.
     *
     * @param kb  on return, \f$ k_{b,j} \f$
     */
    void backwardRateCoefficients(double* const p_kb);
    
    /**
     * Fills the vector ropf with the forward rate of progress variables for
     * each reaction \f$ k_{f,j} \prod_i C_i^{\nu_{ij}^{'}} \Theta_{TB} \f$.
     *
     * @param p_ropf  on return, the forward rates of progress in mol/m^3-s
     */
    void forwardRatesOfProgress(double* const p_ropf);

    /**
     * Fills the vector ropf with the forward rate of progress variables for
     * each reaction \f$ k_{f,j} \prod_i C_i^{\nu_{ij}^{'}} \Theta_{TB} \f$.
     *
     * @param p_conc  the species concentration vector in mol/m^3
     * @param p_ropf  on return, the forward rates of progress in mol/m^3-s
     */
    void forwardRatesOfProgress(
        const double* const p_conc, double* const p_ropf);
    
    /**
     * Fills the vector ropb with the backward rates of progress variables for
     * each reaction \f$ k_{b,j} \prod_i C_i^{\nu_{ij}^{"}} \Theta_{TB} \f$.
     *
     * @param ropb  on return, the backward rates of progress in mol/m^3-s
     */
    void backwardRatesOfProgress(double* const p_ropb);

    /**
     * Fills the vector ropb with the backward rates of progress variables for
     * each reaction \f$ k_{b,j} \prod_i C_i^{\nu_{ij}^{"}} \Theta_{TB} \f$.
     *
     * @param conc  the species concentration vector in mol/m^3
     * @param ropb  on return, the backward rates of progress in mol/m^3-s
     */
    void backwardRatesOfProgress(
        const double* const p_conc, double* const p_ropb);
    
    /**
     * Fills the vector rop with the net rates of progress for each reaction
     * \f$ \left[k_{f,j}\prod_i C_i^{\nu_{ij}^{'}}-k_{b,j}\prod_i 
     * C_i^{\nu_{ij}^"} \right] \Theta_{TB} \f$.
     *
     * @param p_rop on return, the net rates of progress in mol/m^3-s
     */
    void netRatesOfProgress(double* const p_rop);

    /**
     * Fills the vector rop with the net rates of progress for each reaction
     * \f$ \left[k_{f,j}\prod_i C_i^{\nu_{ij}^{'}}-k_{b,j}\prod_i 
     * C_i^{\nu_{ij}^"} \right] \Theta_{TB} \f$.
     *
     * @param p_conc the species concentration vector in mol/m^3
     * @param p_rop on return, the net rates of progress in mol/m^3-s
     */
    void netRatesOfProgress(const double* const p_conc, double* const p_rop);
    
    /**
     * Fills the vector wdot with the net species production rates due to the
     * chemical reactions.
     * \f[
     * \dot{\omega}_i = M_{w,i} \sum_j \left( \nu_{ij}^" - \nu_{ij}^{'} \right) 
     *    \left[ k_{f,j} \prod_i C_i^{\nu_{ij}^{'}} - k_{b,j} \prod_i 
     *    C_i^{\nu_{ij}^"} \right] \Theta_{TB}
     * \f]
     *
     * @param p_wdot - on return, the species production rates in kg/m^3-s
     */
    void netProductionRates(double* const p_wdot);

    /**
     * Fills the matrix p_jac with the species production rate jacobian matrix
     * \f[
     * J_{ij} = \frac{\partial \dot{\omega}_i}{\partial \rho_j}
     * \f]
     * The Jacobian matrix should be sized to be at least the square of the
     * number of species, ns.  Access the jacobian using row-major ordering (ie:
     * J_{ij} = p_jac[i*ns + j]).
     *
     * @param p_jac  - on return, the jacobian matrix \f$J_{ij}\f$
     */
    void jacobianRho(double* const p_jac);

    /**
     * Returns the change in some species quantity across each reaction.
     */
    //void getReactionDelta(
    //    const double* const p_s, double* const p_r) const;


private:

    /**
     * Adds a new reaction to the Kinetics object.
     */
    void addReaction(const Reaction &reaction);
    
    /**
     * Lets the Kinetics object know that the user is done adding reactions
     * allowing the object to optimize how it manages its resources.  Note that
     * calling addReaction() after closeReactions() results in an error.  This
     * method must be called before any other method can be used (except for
     * addReaction()).
     */
    void closeReactions(const bool validate_mechanism = false);

private:

    std::string m_name;

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    
    std::vector<Reaction> m_reactions;

    StoichiometryManager m_reactants;
    StoichiometryManager m_rev_prods;
    StoichiometryManager m_irr_prods;
    
    RateManager*     mp_rates;
    ThirdbodyManager m_thirdbodies;
    JacobianManager  m_jacobian;
    
    double* mp_ropf;
    double* mp_ropb;
    double* mp_rop;
    double* mp_wdot;
};


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_H
