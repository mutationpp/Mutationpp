#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawLH : public GSIRateLaw {

public:
    GSIRateLawLH( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          m_surf_props ( l_data_gsi_rate_law.s_surf_props ) {
    

    }

//=============================================================================================================

    ~GSIRateLawLH( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const { 

       return 0.;

}

//=============================================================================================================

private:
    const SurfaceProperties& m_surf_props;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawLH, GSIRateLaw> gsi_rate_law_langmuir_hinschelwood("L-H");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
