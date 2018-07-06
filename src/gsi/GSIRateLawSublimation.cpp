#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

using namespace Eigen;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawSublimation : public GSIRateLaw
{
public:
    GSIRateLawSublimation(ARGS args)
        : GSIRateLaw(args),
          mv_prod(args.s_products),
          set_state_with_rhoi_T(1),
          pos_T_trans(0),
          idx_gas_prod(0)
    {
        assert(args.s_node_rate_law.tag() == "sublimation");

        args.s_node_rate_law.getAttribute("vap_coef", m_vap_coef,
            "The vapor coefficient for the sublimation reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
            "The pre-exponential coefficient for the sublimation reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute("T", m_activ_en,
            "The activation energy for the sublimation reaction "
            "should be provided.");
    }

//==============================================================================

    ~GSIRateLawSublimation(){}

//==============================================================================

    double forwardReactionRateCoefficient(
        const VectorXd& v_rhoi, const VectorXd& v_Twall) const
    {
    	double Twall = v_Twall(pos_T_trans);

        double sat_vap_p = m_pre_exp*std::exp(-m_activ_en/Twall) ;
        double sat_vap_rho = sat_vap_p * m_thermo.speciesMw(
                                 mv_prod[idx_gas_prod])/(RU*Twall);
                
        m_thermo.setState(
            v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);

        double sp_thermal_speed = m_transport.speciesThermalSpeed(
                                      mv_prod[idx_gas_prod]);

        return (sat_vap_rho - v_rhoi(mv_prod[idx_gas_prod]))*m_vap_coef*
                   sp_thermal_speed / 4.
                   / m_thermo.speciesMw(mv_prod[idx_gas_prod]);
    }

private:
    const size_t pos_T_trans;
    const size_t idx_gas_prod;
    const int set_state_with_rhoi_T;

    double m_vap_coef;
    double m_pre_exp;
    double m_activ_en;

    const std::vector<int>& mv_prod;
};

ObjectProvider<
    GSIRateLawSublimation, GSIRateLaw>
    gsi_rate_law_sublimation("sublimation");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
