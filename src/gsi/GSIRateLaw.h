#ifndef GSI_RATE_LAW_H
#define GSI_RATE_LAW_H

#include "DataGSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Abstract base class that computes the forward reaction rate coefficient
 * for all categories of surface interaction processes.
 */
class GSIRateLaw
{
public:
	/**
	 * Required for self registering different rate laws.
	 */
    typedef const DataGSIRateLaw& ARGS;

    /**
     * Constructor of the abstract base class, taking as input a structure
     * of type DataGSIRateLaw.
     */
    GSIRateLaw(ARGS args)
        : m_thermo(args.s_thermo),
          m_transport(args.s_transport)
    {}

    /**
     * Destructor
     */
    virtual ~GSIRateLaw(){}

    /**
     * Purely virtual function which returns the forward reaction rate
     * coefficients. Its units depend on the specific process and gas-surface
     * interaction model.
     */
    virtual double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const = 0;

protected:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const Mutation::Transport::Transport& m_transport;

}; // class GSIRateLaw

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSI_RATE_LAW_H
