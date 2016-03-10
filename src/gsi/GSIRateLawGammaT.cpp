#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawGammaT : public GSIRateLaw {

public:
    GSIRateLawGammaT( ARGS l_data_gsi_rate_law ) // @todo Extend it to catalysis?
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          v_reactants( l_data_gsi_rate_law.s_reactants ) {

    assert( l_data_gsi_rate_law.s_node_rate_law.tag() == "gamma_T" );

    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
                                      "Error. Nothing provided" );
    l_data_gsi_rate_law.s_node_rate_law.getAttribute( "T", m_activation_energy,
                                      "Error. Nothing provided" );

    }

//=============================================================================================================

    ~GSIRateLawGammaT( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const {
    	double l_Twall = v_Twall(0);

        return computeAverageThermalSpeedforSpeciesI( v_reactants[0], v_Twall ) / 4.0   // @todo change 0;
        		* m_pre_exp * std::exp(- m_activation_energy / l_Twall )
                / m_thermo.speciesMw( v_reactants[0] ) * v_rhoi[0] ;

    }

//=============================================================================================================

private:
    double m_pre_exp;
    double m_activation_energy;

    const std::vector<int>& v_reactants;

//=============================================================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawGammaT, GSIRateLaw> gsi_rate_law_gamma_T("gamma_T");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
