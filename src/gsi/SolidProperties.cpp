/**
 * @file SolidProperties.cpp
 *
 * @brief SolidProperties class when no surface properties are required
 *        to be stored.
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


#include "Utilities.h"

#include "SolidProperties.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

SolidProperties::SolidProperties(
    const Mutation::Utilities::IO::XmlElement& xml_solid_props)
        : m_phi(1.),
          m_h_v(0.)
    {
        if (xml_solid_props.tag() == "solid_properties"){
            xml_solid_props.getAttribute(
                "virgin_to_char_density_ratio", m_phi, 1.);
            xml_solid_props.getAttribute(
                "enthalpy_virgin", m_h_v, 0.);
        }
    }


    } // namespace GasSurfaceInteraction
} // namespace Mutation
