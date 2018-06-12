#include "Thermodynamics.h"
#include "Transport.h"

#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawSublimation : public GSIRateLaw
{
public:
    GSIRateLawSublimation(ARGS args)
        : GSIRateLaw(args),
          v_products(args.s_products)
    {
        assert(args.s_node_rate_law.tag() == "sublimation");

        args.s_node_rate_law.getAttribute( "vap_coef", m_vap_coef,
                                          "Error. Nothing provided" );
        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
                                          "Error. Nothing provided" );
        args.s_node_rate_law.getAttribute( "T", m_activation_energy,
                                          "Error. Nothing provided" );
    }

//=============================================================================================================

    ~GSIRateLawSublimation(){}

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	double l_Twall = v_Twall(0);

        double l_saturated_vapor_p = m_pre_exp * std::exp(-m_activation_energy / l_Twall) ;
        double l_saturated_vapor_rho = l_saturated_vapor_p * m_thermo.speciesMw( v_products[0] ) /
                                       ( Mutation::RU * l_Twall );
                
    	int set_state_with_rhoi_T = 1;
        m_thermo.setState(v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        double l_sp_thermal_speed = m_transport.speciesThermalSpeed(v_products[0]);

        return (l_saturated_vapor_rho - v_rhoi( v_products[0]) ) * m_vap_coef *
               l_sp_thermal_speed / 4
               / m_thermo.speciesMw(v_products[0]);
    }

//=============================================================================================================

private:
    double m_vap_coef;
    double m_pre_exp;
    double m_activation_energy;

    const std::vector<int>& v_products;

//=============================================================================================================

}; // class GSIRateLawSublimation

Mutation::Utilities::Config::ObjectProvider<GSIRateLawSublimation, GSIRateLaw> gsi_rate_law_sublimation("sublimation");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
