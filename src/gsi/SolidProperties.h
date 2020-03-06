/**
 * @file SolidProperties.h
 *
 * @brief  Purely virtual class containing the solid properties.
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


#ifndef SOLID_PROPERTIES_H
#define SOLID_PROPERTIES_H

#include <Eigen/Dense>

namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Structure which stores the necessary inputs for the SolidProperties class.
 */
struct DataSolidProperties
{
    const Mutation::Utilities::IO::XmlElement& s_node_solid_props;
};

//==============================================================================

class SolidProperties
{
public:
    /**
     * Structure containg the information needed by the SolidProperties
     * classes.
     */
    typedef const DataSolidProperties& ARGS;

//==============================================================================

    /**
     * Returns name of this type.
	 */
	static std::string typeName() { return "SolidProperties"; }

//==============================================================================

    /**
     * Default Constructor.
     */
    SolidProperties(ARGS args){}

//==============================================================================

    /**
     * Default Destructor.
     */
    ~SolidProperties(){}

//==============================================================================

    /**
     * Virtual function which returns parameter phi defined as the ratio
     * between the virgin material density and the surface density minus 1.
     */
    virtual double getPhiRatio() const { return 1.; }

//==============================================================================

    /**
     * Virtual function which returns the enthalpy of the virgin material
     * equal to the one set in the gsi input file.
     */
    virtual double getEnthalpyVirginMaterial() const { return 0.; }
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SOLID_PROPERTIES_H
