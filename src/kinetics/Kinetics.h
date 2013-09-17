#ifndef KINETICS_H
#define KINETICS_H

#include "StoichiometryManager.h"
#include "ThirdBodyManager.h"
#include "RateManager.h"
#include "JacobianManager.h"
#include "Reaction.h"
#include "Thermodynamics.h"
#include "Numerics.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Manages the computation of chemical source terms for an entire reaction
 * mechanism.
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
    
    ~Kinetics()
    {
        if (mp_g != NULL) delete [] mp_g;
    }
    
    /**
     * Returns the number of reactions in the mechanism.
     */
    size_t nReactions() const {
        return m_num_rxns;
    }
    
    /**
     * Returns the vector of Reaction objects associated with this kinetics
     * manager.
     */
    const std::vector<Reaction>& reactions() const {
        return m_reactions;
    }
    
    /**
     * Fills the vector keq with the equilibrium constant for each reaction 
     * \f$ K_{C,j} = \left(\frac{P_{atm}}{R_u T}\right)^{\Delta \nu_j}
     * \exp\left(-\frac{\Delta {G^\circ}_j}{R_u T}\right) \f$.
     *
     * @param T    the temperature in K
     * @param keq  on return, \f$ K_{C,j}(T) \f$
     */
    void equilibriumConstants(
        const double T, Mutation::Numerics::RealVector& keq);
    
    /**
     * Fills the vector kf with the forward rate coefficients \f$ k_{f,j} \f$
     * for each reaction.  The coefficients are computed based on the rate laws
     * assigned to each individual reaction.
     *
     * @param T   the temperature in K
     * @param kf  on return, \f$ k_{f,j}(T) \f$
     */
    void forwardRateCoefficients(
        const double T, Mutation::Numerics::RealVector& kf);
    
    /**
     * Fills the vector kb with the backward rate coefficients 
     * \f$ k_{b,j} = k_{f,j} / K_{C,j} \f$ for each reaction.
     *
     * @param T   the temperature in K
     * @param kb  on return, \f$ k_{b,j}(T) \f$
     */
    void backwardRateCoefficients(
        const double T, Mutation::Numerics::RealVector& kb);
    
    /**
     * Fills the vector ropf with the forward rate of progress variables for
     * each reaction \f$ k_{f,j} \prod_i C_i^{\nu_{ij}^{'}} \Theta_{TB} \f$.
     *
     * @param T     the temperature in K
     * @param conc  the species concentration vector in mol/m^3
     * @param ropf  on return, the forward rates of progress in mol/m^3-s
     */
    void forwardRatesOfProgress(
        const double T, const Mutation::Numerics::RealVector& conc, 
        Mutation::Numerics::RealVector& ropf);
    
    /**
     * Fills the vector ropb with the backward rates of progress variables for
     * each reaction \f$ k_{b,j} \prod_i C_i^{\nu_{ij}^{"}} \Theta_{TB} \f$.
     *
     * @param T     the temperature in K
     * @param conc  the species concentration vector in mol/m^3
     * @param ropb  on return, the backward rates of progress in mol/m^3-s
     */
    void backwardRatesOfProgress(
        const double T, const Mutation::Numerics::RealVector& conc, 
        Mutation::Numerics::RealVector& ropb);
    
    /**
     * Fills the vector rop with the net rates of progress for each reaction
     * \f$ \left[k_{f,j}\prod_i C_i^{\nu_{ij}^{'}}-k_{b,j}\prod_i 
     * C_i^{\nu_{ij}^"} \right] \Theta_{TB} \f$.
     *
     * @param T     the temperature in K
     * @param conc  the species concentration vector in mol/m^3
     * @param rop   on return, the net rates of progress in mol/m^3-s
     */
    void netRatesOfProgress(
        const double T, const Mutation::Numerics::RealVector& conc, 
        Mutation::Numerics::RealVector& rop);
    
    /**
     * Fills the vector wdot with the net species production rates due to the
     * chemical reactions.
     * \f[
     * \dot{\omega}_i = M_{w,i} \sum_j \left( \nu_{ij}^" - \nu_{ij}^{'} \right) 
     *    \left[ k_{f,j} \prod_i C_i^{\nu_{ij}^{'}} - k_{b,j} \prod_i 
     *    C_i^{\nu_{ij}^"} \right] \Theta_{TB}
     * \f]
     *
     * @param T     the temperature in K
     * @param conc  the species concentration vector in kg/m^3
     * @param wdot  on return, the species production rates in kg/m^3-s
     */
    void netProductionRates(
        const double T, const double* const conc, double* const wdot);
        
    void netProductionRates(double* const wdot);

    /**
     * Fills the matrix p_jac with the species production rate jacobian matrix
     * \f[
     * J_{ij} = \frac{\partial \dot{\omega}_i}{\partial \rho_j}
     * \f]
     *
     * @param T      - the temperature in K
     * @param p_conc - the species concentration vector in mol/m^3
     * @param p_jac  - on return, the jacobian matrix \f$J_ij\f$
     */
    void jacobianRho(
        const double T, const double *const p_conc, double* const p_jac);

    /**
     * Returns the change in some species quantity across each reaction.
     */
    void getReactionDelta(
        const Mutation::Numerics::RealVector& s, 
        Mutation::Numerics::RealVector& r);
    
    /**
     * @todo Fill in this method.
     */
    double omegaVT();

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

    /**
     * Update quantities that only depend on temperature only if the temperature
     * has changed significantly since the last update.
     */
    void updateT(const double T);

    /**
     * Determins the species index of each species listed in the set of strings.
     */
    //std::vector<size_t> speciesIndices(const std::multiset<std::string>& set);
    
    /**
     * Converts the (species name, efficiency) pairs into (species index,
     * efficiency) pairs.
     */
    //std::vector<std::pair<size_t, double> > thirdbodyEffs(
    //    const std::vector<std::pair<std::string, double> >& string_effs);

private:

    size_t m_num_rxns;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    
    std::vector<Reaction> m_reactions;

    StoichiometryManager m_reactants;
    StoichiometryManager m_rev_prods;
    StoichiometryManager m_irr_prods;
    
    RateManager m_rates;
    ThirdbodyManager m_thirdbodies;
    JacobianManager m_jacobian;
    
    Mutation::Numerics::RealVector m_dnu;
    double* mp_g;
    Mutation::Numerics::RealVector m_lnkf;
    Mutation::Numerics::RealVector m_lnkeq;
    Mutation::Numerics::RealVector m_ropf;
    Mutation::Numerics::RealVector m_ropb;
    Mutation::Numerics::RealVector m_rop;
    
    double m_T_last;
};


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_H
