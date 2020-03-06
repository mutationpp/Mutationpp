/**
 * @file SurfaceRadiation.h
 *
 * @brief Declaration of SurfaceRadiation class.
 */

/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
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


#ifndef SURFACE_RADIATION_H
#define SURFACE_RADIATION_H

#include <Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceState;

/**
 * Class responsible for computing the radiative heat flux of the surface
 * assuming it is grey body. The incoming gas radiation can be considered
 * as an input to the balance by setting the surface feature
 * gas_phase_radiation to true and by calling the gasRadiativeHeatFlux
 * function. If the option gas_phase_radiation is set to false, then
 * the gas phase radiation is either completely neglected or it can
 * be considered equal to the black body radiation at a temperature,
 * an assumption which is not recommended.
 */
class SurfaceRadiation
{
public:
    /**
     *  Constructor
     */
    SurfaceRadiation(
        Mutation::Thermodynamics::Thermodynamics& thermo,
        const Mutation::Utilities::IO::XmlElement& xml_surf_rad,
        const SurfaceState& surf_state,
        bool gas_rad_on);

//==============================================================================
    /**
     *  Destructor
     */
    ~SurfaceRadiation();

//==============================================================================
    /**
     *  Function which returns the net radiative heat flux.
     */
    double surfaceNetRadiativeHeatFlux();

//==============================================================================
    /**
     *  Function which takes the incoming radiative heat from the gas .
     */
    void gasRadiativeHeatFlux(const double& gas_rad_heat_flux){
        if (is_gas_rad_on) m_gas_rad_heat_flux = gas_rad_heat_flux;
    }

//==============================================================================

private:
    const int pos_T_trans;
    const int pos_E;

    const bool is_gas_rad_on;

    double m_eps;
    double m_T_env;
    double m_gas_rad_heat_flux;

    const double m_stef_bolt_const;

    const SurfaceState& m_surf_state;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_RADIATION_H
