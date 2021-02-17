/**
 * @file SolidPropertiesSteadyState.cpp
 *
 * @brief SolidProperties class when solid properties are required
 *        for steady state pyrolysis and solid conduction.
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

#include "AutoRegistration.h"
#include "Utilities.h"

#include "SolidProperties.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SolidPropertiesSteadyState : public SolidProperties
{
public:
    /**
     * Default Constructor
     */
    SolidPropertiesSteadyState(ARGS args)
        : SolidProperties(args),
          m_phi(1.),
          m_h_v(0.)
    {
        if (args.s_node_solid_props.tag() == "solid_properties") {
            args.s_node_solid_props.getAttribute(
                "virgin_to_surf_density_ratio", m_phi, 1.);
            args.s_node_solid_props.getAttribute(
                "enthalpy_virgin", m_h_v, 0.);
        }
    }

//==============================================================================

    /**
     * Default Destructor
     */
    ~SolidPropertiesSteadyState(){}

//==============================================================================
    /**
     * Returns parameter phi defined as the ration between the virgin
     * material density and the surface density minus 1. It is need for the
     * steady state conduction and pyrolysis cases. If it is undefined its
     * value is aytomatically set equal to one.
     */
    double getPhiRatio() const { return 1.; }

//==============================================================================

    /**
     * Returns the enthalpy of the virgin material equal to the one set
     * in the gsi input file. If it is unset this value defaults to 0.
     */
    double getEnthalpyVirginMaterial() const { return 0.; }

private:
    double m_phi;
    double m_h_v;
};

ObjectProvider<
    SolidPropertiesSteadyState, SolidProperties>
    solid_properties_steady_state("steady_state");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
