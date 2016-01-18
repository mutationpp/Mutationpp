#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawAdsorption : public GSIRateLaw {

public:
    GSIRateLawAdsorption( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ){

//    l_data_gsi_rate_law.getAttribute( "sticking_coef", m_S_coef, 
//                                      "The sticking coefficient should be provided for every adsorption reaction. Please provide a sticking_coef!" );
//    l_data_gsi_rate_law.getAttribute( "b", m_beta, 
//                                      "A beta coefficient should be provided for every adsorption reaction. Please provide a 'b' attribute!" );
//    l_data_gsi_rate_law.getAttribute( "E_activation", m_E_act, 
//                                      "Activation energy should be provided for every adsorption reaction. Please provide a 'E_activation' attribute!" );

    

}

//=============================================================================================================

    ~GSIRateLawAdsorption( ){}

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const {

//        double k = m_S_coef * pow( v_Twall(0), m_beta ) * exp( -( m_E_act ) / ( Mutation::RU * v_Twall(0);
//        int idx = 0;

        return 0.0;

//        return k / ( 4 * ( computeAverageThermalSpeedforSpeciesI( idx, v_Twall ) *  total_phase_ mol...
// ) ;
    } 

//=============================================================================================================

private:
    double m_S_coeff;
    double m_beta;
    double m_E_act;

    int m_site;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawAdsorption, GSIRateLaw> gsi_rate_law_adsorption("adsorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
