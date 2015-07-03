#ifndef GASSURFACEINTERACTION_H
#define GASSURFACEINTERACTION_H

#include "Thermodynamics.h"
#include "Transport.h"

#include "SurfaceProperties.h" 
#include "DrivingForces.h" 
#include "SurfaceBalanceSolver.h" 
#include "SurfaceProductionTerms.h"


namespace Mutation {
    namespace GasSurfaceInteraction {

class GasSurfaceInteraction {

public:
    GasSurfaceInteraction( Mutation::Thermodynamics::Thermodynamics& thermo, Mutation::Transport::Transport& transport, const std::string& gsi_catalysis_mechanism_file ) { } 
    ~GasSurfaceInteraction(){ }

    void surfaceProductionRates(){ }
    void solveSurfaceBalance(){ }

private:
//    Mutation::Thermodynamics::Thermodynamics& m_thermo;
//    Mutation::Transport::Transport& m_transport;

    SurfaceProperties m_surf_props;
    DrivingForces m_driving_forces;
    SurfaceBalanceSolver m_surf_solver;
   
    std::vector<SurfaceProductionTerms*> v_surf_production;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GASSURFACEINTERACTION_H
