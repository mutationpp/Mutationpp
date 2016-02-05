#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawAdsorption : public GSIRateLaw {

public:
    GSIRateLawAdsorption( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          m_surf_props ( l_data_gsi_rate_law.s_surf_props ) {

    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "sticking_coef", m_S_coeff, 
                                      "The sticking coefficient should be provided for every adsorption reaction. Please provide a sticking_coef!" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "b", m_beta, 
                                      "A beta coefficient should be provided for every adsorption reaction. Please provide a 'b' attribute!" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "E_activation", m_E_act, 
                                      "Activation energy should be provided for every adsorption reaction. Please provide a 'E_activation' attribute!" );

    const int location_reac_gas_sp = 0;
    m_idx_ads_sp = l_data_gsi_rate_law.s_reactants[location_reac_gas_sp]; 

    }

//=============================================================================================================

    ~GSIRateLawAdsorption( ){}

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const {

    	// @todo Not checked for cases where m_E_act /= 0
        double k = m_S_coeff * pow( v_Twall(0), m_beta ) * exp( -( m_E_act ) / ( Mutation::RU * v_Twall(0) ) );

       // @todo Pay attention when you have different phases!
       // @todo m_surf_props.nTotalSites()^ns!
       return k / ( 4 * computeAverageThermalSpeedforSpeciesI( m_idx_ads_sp, v_Twall ) * m_surf_props.nTotalSites() ) ; // PAY ATTENTION FO R THE UNITS;
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
