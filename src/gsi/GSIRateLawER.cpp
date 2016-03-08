#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawER : public GSIRateLaw {

public:
    GSIRateLawER( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          m_surf_props ( l_data_gsi_rate_law.s_surf_props ) {
    
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "steric_factor", m_S_steric, 
                                      "The steric factor should be provided for every E-R reaction. Please provide a steric factor!" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "b", m_beta, 
                                      "A beta coefficient should be provided for every E-R reaction. Please provide a 'b' attribute!" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "E_activ_recomb", m_E_act, 
                                      "Activation energy should be provided for every E-R reaction. Please provide a 'E_activ_recomb' attribute!" );
    
    const int location_reac_gas_sp = 0;
    m_idx_ads_sp = l_data_gsi_rate_law.s_reactants[location_reac_gas_sp]; 
    

    }

//=============================================================================================================

    ~GSIRateLawER( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const { 

    	int pos_T_trans = 0;
    	double l_T_trans = v_Twall(pos_T_trans);
        double k = m_S_steric * pow( l_T_trans, m_beta ) * exp( -( m_E_act ) / ( Mutation::RU * l_T_trans ) );

       // @todo Pay attention when you have different phases!
       // @todo m_surf_props.nTotalSites()^ns!
       return k / ( 4 * computeAverageThermalSpeedforSpeciesI( m_idx_ads_sp, v_Twall ) * m_surf_props.nTotalSites() ) ; // PAY ATTENTION FOR THE UNITS;

}

//=============================================================================================================

private:
    double m_S_steric;
    double m_beta;
    double m_E_act;

    const SurfaceProperties& m_surf_props;
    int m_idx_ads_sp;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawER, GSIRateLaw> gsi_rate_law_eley_rideal("E-R");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
