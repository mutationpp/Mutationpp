#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawDesorption : public GSIRateLaw {

public:
    GSIRateLawDesorption( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ){ }

//=============================================================================================================

    ~GSIRateLawDesorption( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const { }

//=============================================================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawDesorption, GSIRateLaw> gsi_rate_law_desorption("desorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
