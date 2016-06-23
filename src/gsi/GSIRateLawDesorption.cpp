#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawDesorption : public GSIRateLaw
{
public:
    GSIRateLawDesorption( ARGS args )
        : GSIRateLaw (args)
    {
    	const double l_zero = 0.0;
    	const double l_one = 1.0;

        args.s_node_rate_law.getAttribute("steric_factor", m_steric_factor, l_one);
        args.s_node_rate_law.getAttribute("b", m_beta, l_zero);
        args.s_node_rate_law.getAttribute("freq_vibration_perp", m_vib_perp, l_zero);
        args.s_node_rate_law.getAttribute("E_desorption", m_E_des,
                     "Activation energy should be provided for every desorption reaction. Please provide a 'E_desorption' attribute!");

    }

//=============================================================================================================

    ~GSIRateLawDesorption(){}

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	const double l_low_number = 1.e-30;

    	const size_t pos_T_trans = 0;
    	double l_T_trans = v_Twall(pos_T_trans);

    	double l_freq_factor = m_vib_perp;
    	if (m_vib_perp < l_low_number){
    	    l_freq_factor = Mutation::KB * l_T_trans / Mutation::HP;
    	}

        return (m_steric_factor * m_vib_perp * pow(l_T_trans, m_beta)
               * exp(-(m_E_des)/(Mutation::RU * l_T_trans)));
    }

//=============================================================================================================

private:
    double m_steric_factor;
    double m_beta;
    double m_vib_perp;
    double m_E_des;

}; // class GSIRateLawDesorption

Mutation::Utilities::Config::ObjectProvider<GSIRateLawDesorption, GSIRateLaw> gsi_rate_law_desorption("desorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
