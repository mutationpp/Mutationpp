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
     * @para v_mole_frac_edge
     * @para dx ds
     */
    void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx);

//==============================================================================
/*
 *
 *
 * @param Input mole fractions at the wall
 * @param Output diffusion velocities
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
