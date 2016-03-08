#include "MassBlowingRate.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateNull : public MassBlowingRate {
public:
	MassBlowingRateNull( ARGS l_data_mass_blowing_rate ){ }
	~MassBlowingRateNull(){ }

private:
	double computeBlowingFlux(){ return 0.0; }

};

Mutation::Utilities::Config::ObjectProvider<MassBlowingRateNull, MassBlowingRate> mass_blowing_rate_null("zero"); // @totuesday

    } // namespace GasSurfaceInteraction
} // namespace Mutation
