#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawSublimation : public GSIRateLaw {

public:
    GSIRateLawSublimation( ARGS l_data_gsi_rate_law ) // @todo Extend it to catalysis?
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          v_products( l_data_gsi_rate_law.s_products ) {

    assert( l_data_gsi_rate_law.s_node_rate_law.tag() == "sublimation" );

    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "vap_coef", m_vap_coef,
                                      "Error. Nothing provided" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
                                      "Error. Nothing provided" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "T", m_activation_energy,
                                      "Error. Nothing provided" );

    }

//=============================================================================================================

    ~GSIRateLawSublimation( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const {
    	double l_Twall = v_Twall(0);

        double l_saturated_vapor_p = m_pre_exp * std::exp( -m_activation_energy / l_Twall ) ;
        double l_saturated_vapor_rho = l_saturated_vapor_p * m_thermo.speciesMw( v_products[0] ) /
                                       ( Mutation::RU * l_Twall );
                
        return (l_saturated_vapor_rho - v_rhoi( v_products[0]) ) * m_vap_coef *
               computeAverageThermalSpeedforSpeciesI( v_products[0], v_Twall ) / 4
               / m_thermo.speciesMw( v_products[0] );

    }

//=============================================================================================================

private:
    double m_vap_coef;
    double m_pre_exp;
    double m_activation_energy;

    const std::vector<int>& v_products;

//=============================================================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawSublimation, GSIRateLaw> gsi_rate_law_sublimation("sublimation");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
