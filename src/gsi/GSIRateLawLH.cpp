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

    	double radius_O = 4.8e-11;
    	double radius_N = 5.6e-11;

        return 0.;
    }

//=============================================================================================================

private:
    const SurfaceProperties& m_surf_props;

}; // class GSIRateLawLH

Mutation::Utilities::Config::ObjectProvider<GSIRateLawLH, GSIRateLaw> gsi_rate_law_langmuir_hinschelwood("L-H");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
