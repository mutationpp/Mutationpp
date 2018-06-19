#ifndef WALL_STATE_H
#define WALL_STATE_H

#include <eigen3/Eigen/Dense>

namespace Mutation {
    namespace GasSurfaceInteraction {

class Thermodynamics;

class SurfaceProperties;

class WallState
{
public:
    WallState(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        const SurfaceProperties& surf_props);

    /**
     * Destructor
     */
    ~WallState();
    
    void setWallState(
        const double* const p_mass,
        const double* const p_energy,
        const int state_var = 1);

    void getWallState(
        double* const p_rhoi,
        double* const p_rhoie,
        const int state_var = 1);

    void setWallRhoi(const double* const p_rhoi);
    void setWallT(const double* const p_T);
    void setWallP(const double& p);

    const Eigen::VectorXd& getWallRhoi() const { return mv_rhoi; }
    const Eigen::VectorXd& getWallT() const { return mv_T; }
    double getWallP() const { return m_p; }

    bool isWallStateSet() const{ return m_is_wall_state_set; }

    // Below are FRC Properties!
    void getNdStateGasSurf(Eigen::VectorXd& v_wall_state) const;

private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceProperties& m_surf_props;

    void initializeSurfState();

    const size_t m_ns;
    const size_t m_nT;
    const int m_ns_surf;

    const int m_set_state_rhoi_T;

    Eigen::VectorXd mv_rhoi;
    Eigen::VectorXd mv_T;
    double m_p;

    Eigen::VectorXd mv_surf_props_state;

    bool m_is_wall_state_set;

};

//==============================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // WALL_STATE_H
