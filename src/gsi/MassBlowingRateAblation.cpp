#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "MassBlowingRate.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateAblation : public MassBlowingRate
{
public:
    MassBlowingRateAblation(ARGS args)
        : m_ns(args.s_thermo.nSpecies()),
          v_wall_production_rates(m_ns),
          v_wrk(m_ns+args.s_thermo.nEnergyEqns())
    {
    	std::string prod_term_tag;
        for (size_t i_term = 0; i_term < args.vs_wall_productions_terms.size(); ++i_term){

            prod_term_tag = args.vs_wall_productions_terms[i_term]->getWallProductionTermTag();

            if (prod_term_tag.compare("surface_chemistry") == 0)
                v_wall_production_terms.push_back(args.vs_wall_productions_terms[i_term]);
            if(prod_term_tag.compare("pyrolysis") == 0)
            	v_wall_production_terms.push_back(args.vs_wall_productions_terms[i_term]);
        }

    }

    /**
     * Destructor
     */
    ~MassBlowingRateAblation(){}

    /**
     * This function returns the mass blowing flux in kg/m^2-s as the sum of
     * the heterogeneous reactions.
     */
    double computeBlowingFlux(){

        v_wall_production_rates.setZero();
        v_wrk.setZero();

        for(size_t i_term = 0; i_term < v_wall_production_terms.size(); ++i_term)
        {
        	v_wall_production_terms[i_term]->productionRate(v_wrk);
            v_wall_production_rates += v_wrk.head(m_ns);
        }

        return v_wall_production_rates.sum();
    }

private:
    const size_t m_ns;

    std::vector<WallProductionTerms*> v_wall_production_terms;
    Eigen::VectorXd v_wall_production_rates;
    Eigen::VectorXd v_wrk;

}; //class MassBlowingRateAblation

Mutation::Utilities::Config::ObjectProvider<MassBlowingRateAblation, MassBlowingRate> mass_blowing_rate_ablation("isOn");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
