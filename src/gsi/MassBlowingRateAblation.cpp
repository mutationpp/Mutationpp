#include <eigen3/Eigen/Dense>

#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "MassBlowingRate.h"
#include "WallProductionTerms.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateAblation : public MassBlowingRate
{
public:
    MassBlowingRateAblation(ARGS args)
        : m_ns(args.s_thermo.nSpecies()),
          mv_wall_prod_rates(m_ns),
          mv_wrk(m_ns+args.s_thermo.nEnergyEqns())
    {
    	std::string prod_term_tag;
        for (size_t i_term = 0;
             i_term < args.vs_wall_productions_terms.size();
             ++i_term) {
            prod_term_tag = args.vs_wall_productions_terms[i_term]->
                                getWallProductionTermTag();

            if (prod_term_tag.compare("surface_chemistry") == 0) {
                mv_wall_prod_terms.push_back(
                    args.vs_wall_productions_terms[i_term]);
            }
            if(prod_term_tag.compare("pyrolysis") == 0) {
                mv_wall_prod_terms.push_back(
                    args.vs_wall_productions_terms[i_term]);
            }
        }

    }

//==============================================================================
    /**
     * Destructor
     */
    ~MassBlowingRateAblation(){}

//==============================================================================
    /**
     * This function returns the mass blowing flux in kg/m^2-s as the sum of
     * the heterogeneous reactions.
     */
    double computeBlowingFlux(){

        mv_wall_prod_rates.setZero();
        mv_wrk.setZero();

        for(size_t i_term = 0; i_term < mv_wall_prod_terms.size(); ++i_term)
        {
        	mv_wall_prod_terms[i_term]->productionRate(mv_wrk);
            mv_wall_prod_rates += mv_wrk.head(m_ns);
        }

        return mv_wall_prod_rates.sum();
    }

private:
    const size_t m_ns;

    std::vector<WallProductionTerms*> mv_wall_prod_terms;
    VectorXd mv_wall_prod_rates;
    VectorXd mv_wrk;
};

ObjectProvider<
    MassBlowingRateAblation, MassBlowingRate>
    mass_blowing_rate_ablation("isOn");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
