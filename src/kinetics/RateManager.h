#ifndef KINETICS_RATEMANAGER_H
#define KINETICS_RATEMANAGER_H

#include <vector>

#include "RateLaws.h"
#include "Numerics.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Manages the efficient computation of forward rate coefficients.
 *
 * @see RateLaw
 */
class RateManager
{
public:

    /**
     * Default constructor.
     */
    RateManager() { }
    
    /**
     * Adds a new rate coefficient to the manager.
     *
     * @param rxn   the corresponding reaction index in the Kinetics manager
     * @param rate  the RateLaw to be used
     * 
     * @see RateLaw
     */
    void addRateCoefficient(
        const size_t rxn, const RateLaw *const rate);
    
    /**
     * Returns the forward rate coefficients \f$ k_{f,j}(T) \f$ for each 
     * reaction managed by this RateManager.
     *
     * @param T   temperature in Kelvin
     * @param kf  on return, filled with \f$ k_{f,j}(T) \f$
     */
    void forwardRateCoefficients(
        const double T, Mutation::Numerics::RealVector& kf) const;
    
    /**
     * Returns the natural logarithm of the forward rate coefficients
     * \f$ \ln k_{f,j}(T) \f$.  It is sometimes useful to work with the log
     * of \f$ k_{f,j}(T) \f$ rather than \f$ k_{f,j}(T) \f$ itself.
     *
     * @param T   temperature in Kelvin
     * @param kf  on return, filled with \f$ \ln k_{f,j}(T) \f$
     */
    void lnForwardRateCoefficients(
        const double T, Mutation::Numerics::RealVector& lnkf) const;
    
private:

    std::vector<std::pair<size_t, Arrhenius> > m_arrhenius;

}; // class RateManager


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_RATEMANAGER_H
