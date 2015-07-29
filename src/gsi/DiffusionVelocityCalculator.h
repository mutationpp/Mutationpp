#ifndef DIFFUSIONVELOCITYCALCULATOR_H
#define DIFFUSIONVELOCITYCALCULATOR_H

#include <vector>

#include "Numerics.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "DiffusionDrivingForces.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class DiffusionVelocityCalculator {

public:
    DiffusionVelocityCalculator( const Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport );
    ~DiffusionVelocityCalculator();

    void setDiffusionModel( const Mutation::Numerics::RealVector& lv_mole_frac_edge, const double& l_dx );
    void computeDiffusionVelocities( const Mutation::Numerics::RealVector& lv_mole_frac, Mutation::Numerics::RealVector& lv_diff_velocities );

private:
    Mutation::Transport::Transport& m_transport;

    std::vector<DiffusionDrivingForces*> v_driving_forces;

    Mutation::Numerics::RealVector v_rhs;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSIONVELOCITYCALCULATOR_H
