#include "MassBlowingRate.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateAblation : public MassBlowingRate
{
public:
	MassBlowingRateAblation(ARGS args)
	    : m_wall_production_terms(args.s_wall_productions_terms),
		  v_wall_production_rates(args.s_thermo.nSpecies()){ }

	/**
	 * Destructor
	 */
	~MassBlowingRateAblation(){ }

	/**
	 * This function returns the mass blowing flux in kg/m^2-s as the sum of
	 * the heterogeneous reactions.
	 */
	double computeBlowingFlux(){
		// DO NOT FORGET TO SET THE WALL STATE BEFORE!

		v_wall_production_rates.setZero();
        m_wall_production_terms.productionRate(v_wall_production_rates);

		return v_wall_production_rates.sum();
	}

private:
	WallProductionTerms& m_wall_production_terms; // @todo Check const correctness
    Eigen::VectorXd v_wall_production_rates;

}; //class MassBlowingRateAblation

Mutation::Utilities::Config::ObjectProvider<MassBlowingRateAblation, MassBlowingRate> mass_blowing_rate_ablation("isOn");

    } // namespace GasSurfaceInteraction
} // namespace Mutation

//@todo Blowing Rate Per different solid species.
// REWRITE. Do not include catalytic production rates....
