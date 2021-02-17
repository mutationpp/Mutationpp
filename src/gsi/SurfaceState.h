/**
 * @file SurfaceState.h
 *
 * @brief Declaration of SurfaceState class, responsible for keeping
 * the surface thermodynamic state.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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


#ifndef SURFACE_STATE_H
#define SURFACE_STATE_H

#include <Eigen/Dense>
#include <string>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class DataSolidProperties;
class SurfaceProperties;
class SolidProperties;

//==============================================================================
/*
 * Class which keeps the thermodynamics state of the surface in terms of
 * partial densities and temperatures. The state can be set and retrieved
 * through setter and getter functions. It also acts as an interface for
 * additional surface or bulk properties necessary for gas-surface interaction
 * processes.
 */
class SurfaceState
{
public:
    /**
     * Constructor
     */
    SurfaceState(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        const Mutation::Utilities::IO::XmlElement& xml_surf_props);

    /**
     * Destructor
     */
    ~SurfaceState();

    /**
     * Thermodynamic surface state setter. It sets the state based on the
     * state_var parameter, similarly to how state models work. Currently,
     * setting the state only based on partial densities and temperature is
     * supported (state_var = 1).
     */
    void setSurfaceState(
        const double* const p_mass,
        const double* const p_energy,
        const int state_var = 1);

    /**
     * Thermodynamic surface state getter. It gets the state based on the
     * state_var parameter, similarly to how state models work. Currently,
     * getting the state only based on partial densities and temperature is
     * supported (state_var = 1).
     */
    void getSurfaceState(
        double* const p_rhoi,
        double* const p_rhoie,
        const int state_var = 1) const;

    /**
     * Function which sets the partial densities of the gas phase species
     * on the surface.
     */
    void setSurfaceRhoi(const double* const p_rhoi);

    /**
     * Function which sets the surface temperature.
     */
    void setSurfaceT(const double* const p_T);

    /**
     * Function which gets the partial densities of the gas phase species
     * on the surface.
     */
    const Eigen::VectorXd& getSurfaceRhoi() const { return mv_rhoi; }

    /**
     * Function which gets the surface temperature.
     */
    const Eigen::VectorXd& getSurfaceT() const { return mv_T; }

    /**
     * Returns whether the surface state has been already set.
     */
    bool isSurfaceStateSet() const { return m_is_surface_state_set; }

    /**
     * Helper function to set up the solid properties
     */
    void setSolidProperties(
        const std::string& solid_model,
        const DataSolidProperties& data_solid_props);

    /**
     * Returns a constant reference to the surface properties class.
     */
    const SurfaceProperties& surfaceProps() const {
         return *mp_surf_props;
    }

    /**
     * Returns a constant reference to the surface properties class.
     */
    const SolidProperties& solidProps() const {
        return *mp_solid_props;
    }

private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceProperties* mp_surf_props;
    const SolidProperties* mp_solid_props;

    const size_t m_ns;
    const size_t m_nT;

    Eigen::VectorXd mv_rhoi;
    Eigen::VectorXd mv_T;

    bool m_is_surface_state_set;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_STATE_H
