#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawLH : public GSIRateLaw
{
public:
    GSIRateLawLH(ARGS args)
        : GSIRateLaw(args),
          m_surf_props(args.s_surf_props)
    {}

//=============================================================================================================

    ~GSIRateLawLH(){}

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
        return 0.;
    }

//=============================================================================================================

private:
    const SurfaceProperties& m_surf_props;

}; // class GSIRateLawLH

Mutation::Utilities::Config::ObjectProvider<GSIRateLawLH, GSIRateLaw> gsi_rate_law_langmuir_hinschelwood("L-H");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
