#include "Thermodynamics.h"
#include "Transport.h"

#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawGammaT : public GSIRateLaw
{
public:
    GSIRateLawGammaT(ARGS args)
                        : GSIRateLaw(args),
                          v_reactants(args.s_reactants)
    {
        assert(args.s_node_rate_law.tag() == "gamma_T");

        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
                                          "Error. Nothing provided" );
        args.s_node_rate_law.getAttribute( "T", m_activation_energy,
                                          "Error. Nothing provided" );

    }

//=============================================================================================================

    ~GSIRateLawGammaT( ){ }

//=============================================================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	double l_Twall = v_Twall(0);

    	const int set_state_with_rhoi_T = 1;
    	m_thermo.setState(v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
    	double m_sp_thermal_speed = m_transport.speciesThermalSpeed(v_reactants[0]);  // @todo change 0;

        return  m_sp_thermal_speed/4.0
                * m_pre_exp * std::exp(- m_activation_energy/l_Twall)
                / m_thermo.speciesMw(v_reactants[0])*v_rhoi(v_reactants[0]);
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
