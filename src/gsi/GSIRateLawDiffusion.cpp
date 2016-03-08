#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawDiffusion : public GSIRateLaw {

public:
    GSIRateLawDiffusion( ARGS l_data_gsi_rate_law )
                        : GSIRateLaw ( l_data_gsi_rate_law ),
                          m_surf_props ( l_data_gsi_rate_law.s_surf_props ) {
    

    }

//=============================================================================================================

    ~GSIRateLawDiffusion( ){ }

//=============================================================================================================

    double forwardReactionRate( const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall ) const { 

       return 0.;

}

//=============================================================================================================

private:
    const SurfaceProperties& m_surf_props;

};

Mutation::Utilities::Config::ObjectProvider<GSIRateLawDiffusion, GSIRateLaw> gsi_rate_law_diffusion("SurfaceDiffusion");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
