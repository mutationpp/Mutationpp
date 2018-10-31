/**
 * @file SolidProperties.h
 *
 * @brief  Purely virtual class SurfaceProperties.
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


#ifndef SOLID_PROPERTIES_H
#define SOLID_PROPERTIES_H

#include <eigen3/Eigen/Dense>


namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SolidProperties
{
public:
    /**
     * Default Constructor.
     */
    SolidProperties(
        const Mutation::Utilities::IO::XmlElement& xml_solid_props);

//==============================================================================

    /**
     * Default Destructor.
     */
    ~SolidProperties(){}

//==============================================================================

    /**
     *
     */
    double getPhiRatio() const { return m_phi; }
//==============================================================================

    /**
     *
     */
    double getEnthalpyVirginMaterial() const { return m_h_v; }

//==============================================================================
private:
    double m_phi;
    double m_h_v;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SOLID_PROPERTIES_H
