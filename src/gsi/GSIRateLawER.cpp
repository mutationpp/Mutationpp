#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawER : public GSIRateLaw
{
public:
    GSIRateLawER(ARGS args)
        : GSIRateLaw (args),
          m_surf_props (args.s_surf_props)
    {
        const double l_huge_number = 1.e20;
        const double l_zero = 0.0;

        args.s_node_rate_law.getAttribute( "steric_factor", m_S_steric,
                                      "The steric factor should be provided for every E-R reaction. Please provide a steric factor!" );
        args.s_node_rate_law.getAttribute( "b", m_beta,
                                      "A beta coefficient should be provided for every E-R reaction. Please provide a 'b' attribute!" );
        args.s_node_rate_law.getAttribute( "E_activ_recomb", m_E_act,
                                      "Activation energy should be provided for every E-R reaction. Please provide a 'E_activ_recomb' attribute!" );
    
        args.s_node_rate_law.getAttribute("T_threshold", m_T_threshold, l_huge_number);
        args.s_node_rate_law.getAttribute("sticking_coef_correction", m_beta_T_corr, l_zero);

       const int location_reac_gas_sp = 0;
       m_idx_er_sp = args.s_reactants[location_reac_gas_sp];

    }

//=============================================================================================================

    ~GSIRateLawER(){}

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
        int set_state_with_rhoi_T = 1;
    	const size_t pos_T_trans = 0;
    	double l_T_trans = v_Twall(pos_T_trans);

        m_thermo.setState(v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        double l_sp_thermal_speed = m_transport.speciesThermalSpeed(m_idx_er_sp);

        //  Sticking Coefficient Temperature Correction
        double l_S_steric_corr = m_S_steric;
        if ( m_T_threshold < l_T_trans ) {
        	l_S_steric_corr *= exp( -m_beta_T_corr *(l_T_trans - m_T_threshold) );
        }

        double k = l_S_steric_corr * pow(l_T_trans, m_beta) * exp(-(m_E_act) /
        		   (Mutation::RU * l_T_trans));

       // @todo Pay attention when you have different phases!
       return k / (4 * l_sp_thermal_speed * m_surf_props.nTotalSites()) ;
    }

//=============================================================================================================

private:
    double m_S_steric;
    double m_beta;
    double m_E_act;

    double m_T_threshold;
    double m_beta_T_corr;

    const SurfaceProperties& m_surf_props;
    int m_idx_er_sp;

}; // class GSIRateLawER

Mutation::Utilities::Config::ObjectProvider<GSIRateLawER, GSIRateLaw> gsi_rate_law_eley_rideal("E-R");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
