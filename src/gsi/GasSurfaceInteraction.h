#ifndef GASSURFACEINTERACTION_H
#define GASSURFACEINTERACTION_H

#include <Eigen/Dense>
#include <string>

#include "Thermodynamics.h"
#include "Transport.h"

#include "SurfaceProperties.h"
#include "WallState.h"
#include "SurfaceBalanceSolver.h" 

namespace Mutation {
    namespace GasSurfaceInteraction {

class GasSurfaceInteraction {

public:
    GasSurfaceInteraction( Mutation::Thermodynamics::Thermodynamics& l_thermo, 
                           Mutation::Transport::Transport& l_transport, 
                           std::string l_gsi_mechanism_file );
    ~GasSurfaceInteraction();

    void setWallState( const double* const l_mass, const double* const l_energy, const int state_variable );
    void getWallState( double* const l_mass, double* const l_energy, const int state_variable );

    void surfaceProductionRates( double * const lp_mass_prod_rate );
    
    void setDiffusionModel( const double* const lp_mole_frac_edge, const double& l_dx );
    void solveSurfaceBalance();

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Transport::Transport& m_transport;

    SurfaceProperties* mp_surf_props;
    WallState* mp_wall_state;
    SurfaceBalanceSolver* mp_surf_solver;

    std::string m_gsi_mechanism;

    Eigen::VectorXd v_mass_prod_rate; // Ugly! To fix

    inline void locateGSIInputFile( std::string& l_gsi_mechanism_file );
    inline void errorWrongTypeofGSIFile( const std::string& l_gsi_root_tag );
    inline void errorInvalidGSIFileProperties( const std::string& l_gsi_option );

    // TEMPORARY
    Eigen::VectorXd v_mole_frac_edge;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GASSURFACEINTERACTION_H
