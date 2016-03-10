#include "MassBlowingRate.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateAblation : public MassBlowingRate {
public:
	MassBlowingRateAblation( ARGS l_data_mass_blowing_rate )
	                       : m_wall_production_terms( l_data_mass_blowing_rate.s_wall_productions_terms ),
							 v_wall_production_rates( l_data_mass_blowing_rate.s_thermo.nSpecies() ){ }

	~MassBlowingRateAblation(){ }

	double computeBlowingFlux(){
		// DO NOT FORGET TO SET THE WALL STATE BEFORE!

		v_wall_production_rates.setZero();
        m_wall_production_terms.productionRate( v_wall_production_rates );

		return v_wall_production_rates.sum();
	}

private:
	WallProductionTerms& m_wall_production_terms; // @todo Check const correctness

    Eigen::VectorXd v_wall_production_rates;

};

Mutation::Utilities::Config::ObjectProvider<MassBlowingRateAblation, MassBlowingRate> mass_blowing_rate_ablation("isOn"); // @totuesday

    } // namespace GasSurfaceInteraction
} // namespace Mutation

//@todo Blowing Rate Per different solid species.
// REWRITE. Do not include catalytic production rates....
