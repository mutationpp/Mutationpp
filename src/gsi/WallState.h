/**
 * @file WallState.h
 *
 * @brief Declaration of WallState class.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


#ifndef WALL_STATE_H
#define WALL_STATE_H

#include <eigen3/Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProperties;

//==============================================================================
/*
 *
 */
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

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // WALL_STATE_H
