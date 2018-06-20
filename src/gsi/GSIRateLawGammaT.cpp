#include "Thermodynamics.h"
#include "Transport.h"

#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawGammaT : public GSIRateLaw
{
public:
    GSIRateLawGammaT(ARGS args)
        : GSIRateLaw(args),
          mv_react(args.s_reactants),
          pos_T_trans(0),
          idx_react(0)
    {
        assert(args.s_node_rate_law.tag() == "gamma_T");

        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
            "The pre-exponential coefficient for the reaction "
            "should be provided with gamma as a function of temperature.");
        args.s_node_rate_law.getAttribute( "T", m_activ_en,
            "The activation energy for the reaction "
            "should be provided with gamma as a function of temperature.");
    }

//==============================================================================

    ~GSIRateLawGammaT( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	double Twall = v_Twall(pos_T_trans);

    	const int set_state_with_rhoi_T = 1;
    	m_thermo.setState(
    	    v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
    	double m_sp_thermal_speed = m_transport.speciesThermalSpeed(
                                        mv_react[idx_react]);

        return  m_sp_thermal_speed/4.
                * m_pre_exp * std::exp(- m_activ_en/Twall)
                / m_thermo.speciesMw(
                mv_react[idx_react])*v_rhoi(mv_react[idx_react]);
    }

private:
    const size_t pos_T_trans;
    const size_t idx_react;

    double m_pre_exp;
    double m_activ_en;

    const std::vector<int>& mv_react;
};

ObjectProvider<
    GSIRateLawGammaT, GSIRateLaw>
    gsi_rate_law_gamma_T("gamma_T");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
