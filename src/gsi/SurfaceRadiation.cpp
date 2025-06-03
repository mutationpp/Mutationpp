/**
 * @file SurfaceRadiation.cpp
 *
 * @brief Implementation of the SurfaceRadiation class.
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

#include "Constants.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceRadiation.h"
#include "SurfaceState.h"

using namespace Eigen;

using namespace Mutation;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace GasSurfaceInteraction {

SurfaceRadiation::SurfaceRadiation(
    Mutation::Thermodynamics::Thermodynamics& thermo,
    const Mutation::Utilities::IO::XmlElement& xml_surf_rad,
    const SurfaceState& surf_state,
    bool gas_rad_on)
        :  m_surf_state(surf_state),
           m_gas_rad_heat_flux(0.),
           is_gas_rad_on(gas_rad_on),
           pos_E(thermo.nSpecies()),
           pos_T_trans(0)
{
    xml_surf_rad.getAttribute("emissivity", m_eps,
        "Error in the surface radiation input. Surface emissivity "
        "coefficient should be provided");
    xml_surf_rad.getAttribute("T_env", m_T_env, 0.);

    if (!is_gas_rad_on)
        m_gas_rad_heat_flux = Mutation::SB * std::pow(m_T_env, 4.0);
}

//==============================================================================

SurfaceRadiation::~SurfaceRadiation(){}

//==============================================================================

double SurfaceRadiation::surfaceNetRadiativeHeatFlux()
{
    double T_surf = (m_surf_state.getSurfaceT())(pos_T_trans);
    return m_eps * (Mutation::SB * std::pow(T_surf, 4.0) - m_gas_rad_heat_flux);
}

    } // namespace GasSurfaceInteraction
} // namespace Mutation
