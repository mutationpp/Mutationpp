#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawDiffusion : public GSIRateLaw
{
public:
    GSIRateLawDiffusion(ARGS args)
        : GSIRateLaw (args),
          m_surf_props (args.s_surf_props)
    {}

//=============================================================================================================

    ~GSIRateLawDiffusion(){}

//=============================================================================================================

    double forwardReactionRateCoeffient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
       return 0.;
    }

//=============================================================================================================

private:
    const SurfaceProperties& m_surf_props;

}; // class GSIRateLawDiffusion

Mutation::Utilities::Config::ObjectProvider<GSIRateLawDiffusion, GSIRateLaw> gsi_rate_law_diffusion("SurfaceDiffusion");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
