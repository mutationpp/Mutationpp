#ifndef GSIRATEMANAGER_H
#define GSIRATEMANAGER_H

#include<Eigen/Dense>

#include "DataGSIRateManager.h"
#include "GSIReaction.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateManager {

public:
    typedef DataGSIRateManager ARGS;

    GSIRateManager( ARGS l_data_rate_manager ) 
                  : m_thermo( l_data_rate_manager.s_thermo ),
                    m_surf_props( l_data_rate_manager.s_surf_props ),
                    m_wall_state( l_data_rate_manager.s_wall_state ),
                    v_reactions( l_data_rate_manager.s_reactions ) { }
    virtual ~GSIRateManager(){ };

    virtual void computeRate( Eigen::VectorXd& lv_mass_prod_rate ) = 0;

protected:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceProperties& m_surf_props;
    const WallState& m_wall_state;
    const std::vector<GSIReaction*>& v_reactions;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSIRATEMANAGER_H
