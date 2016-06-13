#ifndef GSI_RATE_MANAGER_H
#define GSI_RATE_MANAGER_H

#include<Eigen/Dense>

#include "DataGSIRateManager.h"
#include "GSIReaction.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateManager
{
public:
	/**
	 * Required for self registering different rate laws.
	 */
    typedef DataGSIRateManager ARGS;

    /**
     * Constructor of the abstract base class, taking as input a structure
     * of type DataGSIRateManager.
     */
    GSIRateManager(ARGS l_data_rate_manager)
       : m_thermo( l_data_rate_manager.s_thermo ),
         m_surf_props( l_data_rate_manager.s_surf_props ),
         m_wall_state( l_data_rate_manager.s_wall_state ),
         v_reactions( l_data_rate_manager.s_reactions ) { }

    /**
     * Destructor
     */
    virtual ~GSIRateManager(){ };

    /**
     * This purely virtual function returns the surface production rates according
     * to the selected gas surface interaction model.
     *
     * @return Vector of the gas phase species production rates in kg/m^2-s
     */
    virtual Eigen::VectorXd computeRate() = 0;

protected:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceProperties& m_surf_props;
    const WallState& m_wall_state;
    const std::vector<GSIReaction*>& v_reactions;

}; // class GSIRateManager

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSI_RATE_MANAGER_H
