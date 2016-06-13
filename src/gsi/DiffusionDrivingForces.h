#ifndef DIFFUSION_DRIVING_FORCES_H
#define DIFFUSION_DRIVING_FORCES_H

#include <Eigen/Dense>

namespace Mutation {
    namespace GasSurfaceInteraction {

// JB: I made several formatting changes to give an example for spacing and
// indention
class DiffusionDrivingForces
{
public:

    DiffusionDrivingForces(
        const Mutation::Thermodynamics::Thermodynamics& l_thermo)
        : v_mole_frac_edge(l_thermo.nSpecies()),
          // JB: don't use 0 and 1 for bools!!
          m_diff_model_set(false)
    { }

    // JB: should comment even empty functions
    ~DiffusionDrivingForces(){ }

    // JB: please format functions with opening bracket on its own line, unless the function
    // is very short (1-2 lines), no extra whitespace lines at the end, and start
    // long argument lists on a new line with a tab
    void computeDrivingForces(
        const Eigen::VectorXd& lv_mole_frac, Eigen::VectorXd& lv_driving_force)
    {
        // JB: if states with more than one code line should be on multiple lines
        if ( m_diff_model_set == 0 ) {
            std::cerr << "Error diffusion model!" << std::endl;
            exit(1);
        }

        lv_driving_force = (lv_mole_frac - v_mole_frac_edge) / m_dx;
    }

    void setDiffusionCalculator(
        const Eigen::VectorXd& lv_mole_frac_edge, const double& l_dx)
    {
        v_mole_frac_edge = lv_mole_frac_edge;
        m_dx = l_dx;
        // JB: don't use 0 and 1 for bools!!
        m_diff_model_set = true;
    }

private:

    Eigen::VectorXd v_mole_frac_edge;
    double m_dx;
    bool m_diff_model_set;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSION_DRIVING_FORCES_H
