#include "AutoRegistration.h"
#include "Thermodynamics.h"

#include "MassBlowingRate.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateNull : public MassBlowingRate
{
public:
	MassBlowingRateNull(ARGS args){ }
	~MassBlowingRateNull(){ }

private:
	/**
	 * Returns blowing flux equal to zero.
	 */
	double computeBlowingFlux(){ return 0.0; }

}; // class MassBlowingRateNull

ObjectProvider<
    MassBlowingRateNull, MassBlowingRate>
    mass_blowing_rate_null("zero");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
