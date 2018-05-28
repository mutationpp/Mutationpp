#ifndef DIFFUSION_DRIVING_FORCES_H
#define DIFFUSION_DRIVING_FORCES_H

#include <eigen3/Eigen/Dense>

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
        m_diff_model_set = true;
    }

    // @TO REMOVE
    void setTemperatureGradientCalculator(
        const double& l_T, const double& l_dx_T)
    {
        m_T_edge = l_T;
        m_dx_T = l_dx_T;
    }

    double computeTemperatureDrivingForce(const double& l_T)
    {
        return -(l_T - m_T_edge) / m_dx_T;
    }

private:

    Eigen::VectorXd v_mole_frac_edge;
    double m_dx;
    bool m_diff_model_set;

    double m_T_edge;
    double m_dx_T;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DIFFUSION_DRIVING_FORCES_H
