#ifndef DIFFUSION_VELOCITY_CALCULATOR_H
#define DIFFUSION_VELOCITY_CALCULATOR_H

namespace Mutation {
    namespace GasSurfaceInteraction {

class VectorXd;

class Thermodynamics;
class Transport;

class DiffusionVelocityCalculator
{
public:
    DiffusionVelocityCalculator(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        Mutation::Transport::Transport& transport);

    ~DiffusionVelocityCalculator();

    void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx);

    void computeDiffusionVelocities(
        const Eigen::VectorXd& v_mole_frac,
        Eigen::VectorXd& v_diff_velocities);

private:
    Mutation::Transport::Transport& m_transport;

    Eigen::VectorXd mv_mole_frac_edge;
    Eigen::VectorXd mv_dxidx;

    double m_dx;
    bool m_is_diff_set;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSION_VELOCITY_CALCULATOR_H
