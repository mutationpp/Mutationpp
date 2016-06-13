#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawAdsorption : public GSIRateLaw
{
public:
    GSIRateLawAdsorption(ARGS args)
        : GSIRateLaw(args),
          m_surf_props (args.s_surf_props)
    {
    args.s_node_rate_law.getAttribute("sticking_coef", m_S_coeff, errorStickingCoef().c_str());
    args.s_node_rate_law.getAttribute("b", m_beta, errorBeta().c_str());
    args.s_node_rate_law.getAttribute("E_activation", m_E_act, errorActivationEnergy().c_str());

    const int location_reac_gas_sp = 0;
    m_idx_ads_sp = args.s_reactants[location_reac_gas_sp];
    }

//=============================================================================================================

    ~GSIRateLawAdsorption(){}

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	int set_state_with_rhoi_T = 1;
        m_thermo.setState(v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        double l_sp_thermal_speed = m_transport.speciesThermalSpeed(m_idx_ads_sp);

    	// @todo Not checked for cases where m_E_act /= 0
        double k = m_S_coeff * pow(v_Twall(0), m_beta) * exp(-(m_E_act) /
        		   (Mutation::RU*v_Twall(0)));

       // @todo Pay attention when you have different phases!
       // @todo m_surf_props.nTotalSites()^ns!

       return k / (4 * l_sp_thermal_speed * m_surf_props.nTotalSites());
    } 

//=============================================================================================================

private:
    const std::string errorStickingCoef() const
    {
        return "The sticking coefficient should be provided for every adsorption reaction. Please provide a sticking_coef!";
    }

    const std::string errorBeta() const
    {
        return "A beta coefficient should be provided for every adsorption reaction. Please provide a 'b' attribute!";
    }

    const std::string errorActivationEnergy() const
    {
    	return "Activation energy should be provided for every adsorption reaction. Please provide a 'E_activation' attribute!";
    }

//=============================================================================================================
private:
    double m_S_coeff;
    double m_beta;
    double m_E_act;

    int m_idx_ads_sp;

    const SurfaceProperties& m_surf_props;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawAdsorption, GSIRateLaw> gsi_rate_law_adsorption("adsorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
