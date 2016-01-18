#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawER : public GSIRateLaw {

public:
    GSIRateLawER( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ) {}

//=============================================================================================================

    ~GSIRateLawER( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const { }

//=============================================================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawER, GSIRateLaw> gsi_rate_law_eley_rideal("E-R");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
