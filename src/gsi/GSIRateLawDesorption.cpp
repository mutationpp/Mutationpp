#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawDesorption : public GSIRateLaw {

public:
    GSIRateLawDesorption( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ){

    l_data_gsi_rate_law.s_node_rate_law.getAttribute("steric_factor", m_steric_factor,
                      "The sticking coefficient should be provided for every desorption reaction. Please provide a sticking_coef!" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute("b", m_beta, 
                      "A beta coefficient should be provided for every desorption reaction. Please provide a 'b' attribute!" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute("freq_vibration_perp", m_vib_perp, "In a thermal desorption reaction the vibrational frequency perpendicular to the surface should be provide");
    l_data_gsi_rate_law.s_node_rate_law.getAttribute("E_desorption", m_E_des, 
                      "Activation energy should be provided for every desorption reaction. Please provide a 'E_desorption' attribute!" );

}

//=============================================================================================================

    ~GSIRateLawDesorption( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const { 
        return( m_steric_factor * m_vib_perp * pow(v_Twall[0], m_beta)
                                  * exp(-(m_E_des)/(RU * v_Twall[0] )) );
    }


//=============================================================================================================

private:
    double m_steric_factor;
    double m_beta;
    double m_vib_perp;
    double m_E_des;


};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawDesorption, GSIRateLaw> gsi_rate_law_desorption("desorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
