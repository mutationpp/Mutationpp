#ifndef DIFFUSION_VELOCITY_CALCULATOR_H
#define DIFFUSION_VELOCITY_CALCULATOR_H

#include <eigen3/Eigen/Dense>

namespace Mutation { namespace Thermodynamics {class Thermodynamics; }}
namespace Mutation { namespace Transport {class Transport; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

/*
 * Class responsible for computing the
 */
class DiffusionVelocityCalculator
{
public:
	/*
	 * Constructor of the class. It takes as arguments a copy of the
	 * thermodynamics class and a copy of the transport class in order
	 * to call the Stefan-Maxwell method.
	 */
    DiffusionVelocityCalculator(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        Mutation::Transport::Transport& transport);

//==============================================================================
    /*
     * Default Destructor
     */
    ~DiffusionVelocityCalculator();

//==============================================================================
    /*
     * Function used to set up the mole fractions at a distance from the wall.
     *
     * @para mole_frac_edge mole fractions of the species at a distance
     * @para dx distance in m
     */
    void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx);

//==============================================================================
    /*
     * Function used to compute the diffusion velocities with the given
     * mole fraction at the wall the ones imposed at a given distance from the
     * wall after the diffusion model has been set. The mole fraction gradients
     * are computed with simple first order differentiation.
     *
     * @param mole_frac mole fractions of the species at the wall
     * @param on return diffusion velocities in m/s
     *
     */
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
