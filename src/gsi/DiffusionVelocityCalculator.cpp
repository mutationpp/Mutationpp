#include "Errors.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "DiffusionVelocityCalculator.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//==============================================================================

DiffusionVelocityCalculator::DiffusionVelocityCalculator(
    const Mutation::Thermodynamics::Thermodynamics& thermo,
    Mutation::Transport::Transport& transport)
        : mv_dxidx(thermo.nSpecies()),
          mv_mole_frac_edge(thermo.nSpecies()),
          m_transport(transport),
          m_dx(0.),
          m_is_diff_set(0)
{ }

//==============================================================================

DiffusionVelocityCalculator::~DiffusionVelocityCalculator(){}

//==============================================================================

void DiffusionVelocityCalculator::setDiffusionModel(
    const Eigen::VectorXd& v_mole_frac_edge, const double& dx)
{
    mv_mole_frac_edge = v_mole_frac_edge;

    if (m_dx <= 0.) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::setDiffusionModel() with a "
        << "distance less or equal to zero. The distance dx should always be "
           "positive";
    }
    m_dx = dx;

    m_is_diff_set = true;
}

//==============================================================================

void DiffusionVelocityCalculator::computeDiffusionVelocities(
    const Eigen::VectorXd& v_mole_frac,
    Eigen::VectorXd& v_diff_velocities)
{
    if (!m_is_diff_set) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::computeDrivingForces() before "
        << "calling DiffusionVelocityCalculator::setDiffusionCalculator()";
    }

    mv_dxidx = (v_mole_frac - mv_mole_frac_edge)/m_dx;

    double electric_field = 0.E0;
    m_transport.stefanMaxwell(
        mv_dxidx.data(),
        v_diff_velocities.data(),
        electric_field);
}

//==============================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
