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
    const double l_huge_number = 1.e20;
    const double l_zero = 0.0;

    args.s_node_rate_law.getAttribute("sticking_coef", m_S_coeff, errorStickingCoef().c_str());
    args.s_node_rate_law.getAttribute("b", m_beta, errorBeta().c_str());
    args.s_node_rate_law.getAttribute("E_activation", m_E_act, errorActivationEnergy().c_str());
    args.s_node_rate_law.getAttribute("T_threshold", m_T_threshold, l_huge_number);
    args.s_node_rate_law.getAttribute("sticking_coef_correction", m_beta_T_corr, l_zero);

    const int location_reac_gas_sp = 0;
    const int n_sticking_species = 1;
    m_idx_ads_sp = args.s_reactants[location_reac_gas_sp];
    m_stick_coef_power = args.s_reactants.size() - n_sticking_species;
    }

//=============================================================================================================

    ~GSIRateLawAdsorption(){}

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	int set_state_with_rhoi_T = 1;
    	const size_t pos_T_trans = 0;
    	double l_T_trans = v_Twall(pos_T_trans);

        m_thermo.setState(v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        double l_sp_thermal_speed = m_transport.speciesThermalSpeed(m_idx_ads_sp);

        //  Sticking Coefficient Temperature Correction
        double l_S_coef_T_corr = m_S_coeff;
        if (m_T_threshold < l_T_trans) {
        	l_S_coef_T_corr *= exp(-m_beta_T_corr * (l_T_trans - m_T_threshold));
        }

    	// @todo Not checked for cases where m_E_act /= 0
        double k = l_S_coef_T_corr * pow(l_T_trans, m_beta) * exp(-(m_E_act) /
        		   (Mutation::RU * l_T_trans));

       // @todo Pay attention when you have different phases!
       return k / (4 * l_sp_thermal_speed * pow(m_surf_props.nTotalSites(), m_stick_coef_power));
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

    double m_T_threshold;
    double m_beta_T_corr;

    int m_idx_ads_sp;
    int m_stick_coef_power;

    const SurfaceProperties& m_surf_props;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawAdsorption, GSIRateLaw> gsi_rate_law_adsorption("adsorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
