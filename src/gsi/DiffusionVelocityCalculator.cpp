#include <stddef.h>

#include "DiffusionVelocityCalculator.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

DiffusionVelocityCalculator::DiffusionVelocityCalculator( const Mutation::Thermodynamics::Thermodynamics& l_thermo, Mutation::Transport::Transport& l_transport )
                                                        : v_rhs( l_thermo.nSpecies() ),
                                                          m_transport( l_transport ) {

    v_driving_forces.push_back( new DiffusionDrivingForces( l_thermo ) );

}

//======================================================================================

DiffusionVelocityCalculator::~DiffusionVelocityCalculator(){
    delete v_driving_forces[0];
}

//======================================================================================

void DiffusionVelocityCalculator::setDiffusionModel(
    const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx)
{
    v_driving_forces[0]->setDiffusionCalculator(lv_mole_frac_edge, l_dx);
}

//======================================================================================

void DiffusionVelocityCalculator::setConductiveHeatFluxModel(
    const double& l_T, const double& l_dx_T)
{
    v_driving_forces[0]->setTemperatureGradientCalculator(l_T, l_dx_T);
}

//======================================================================================

void DiffusionVelocityCalculator::computeDiffusionVelocities(
    const Eigen::VectorXd& lv_mole_frac, Eigen::VectorXd& lv_diff_velocities)
{
    v_rhs.setZero();

    for (int i_driv = 0; i_driv < v_driving_forces.size(); ++i_driv)
    {
        v_driving_forces[i_driv]->computeDrivingForces(lv_mole_frac, v_rhs);
    }

    double electric_field = 0.E0;
    m_transport.stefanMaxwell(v_rhs.data(), lv_diff_velocities.data(), electric_field);
}

//======================================================================================
double DiffusionVelocityCalculator::computeConductiveHeatFlux(const double& l_T)
{
    // The state is already set
    double lambda = m_transport.frozenThermalConductivity();
    return lambda * v_driving_forces[0]->computeTemperatureDrivingForce(l_T);

}
//======================================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
