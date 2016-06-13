#ifndef DIFFUSION_VELOCITY_CALCULATOR_H
#define DIFFUSION_VELOCITY_CALCULATOR_H

#include <vector>
#include <Eigen/Dense>

#include "Thermodynamics.h"
#include "Transport.h"

#include "DiffusionDrivingForces.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class DiffusionVelocityCalculator {

public:
    DiffusionVelocityCalculator( const Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport );
    ~DiffusionVelocityCalculator();

    void setDiffusionModel( const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx );
    void computeDiffusionVelocities( const Eigen::VectorXd& lv_mole_frac, Eigen::VectorXd& lv_diff_velocities );

private:
    Mutation::Transport::Transport& m_transport;

    std::vector<DiffusionDrivingForces*> v_driving_forces;

    Eigen::VectorXd v_rhs; // Rename to something that I will understand

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSION_VELOCITY_CALCULATOR_H
