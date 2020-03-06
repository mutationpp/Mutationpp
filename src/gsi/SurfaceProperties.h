/**
 * @file SurfaceProperties.h
 *
 * @brief  Purely virtual class SurfaceProperties.
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


#ifndef SURFACE_PROPERTIES_H
#define SURFACE_PROPERTIES_H

#include <Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Structure which stores the necessary inputs for the SurfaceProperties class.
 */
struct DataSurfaceProperties
{
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Utilities::IO::XmlElement& s_node_surf_props;
};

//==============================================================================

class SurfaceProperties
{
public:
	/**
	 * Structure containg the information needed by the SurfaceProperties
     * classes.
	 */
    typedef const DataSurfaceProperties& ARGS;

//==============================================================================
    /**
     * Returns name of this type.
	 */
	static std::string typeName() { return "SurfaceProperties"; }

//==============================================================================

    /**
     * Default Constructor.
     */
    SurfaceProperties(ARGS args){ }

//==============================================================================

    /**
     * Default Destructor.
     */
    virtual ~SurfaceProperties(){ }

//==============================================================================
    /**
     * Returns the index associated with the surface species. It is always
     * number of gas species plus the surface phase index.
     */
    virtual int surfaceSpeciesIndex(const std::string& str_sp) const {
        return -1;
    }

//==============================================================================
    /**
     * Returns the gas phase species associated with the surface species.
     */
    virtual int surfaceToGasIndex(const int& i_surf_species) const {
        return -1;
    }

//==============================================================================
    /**
     * Returns the number of surface species.
     */
    virtual size_t nSurfaceSpecies() const { return 0; }

//==============================================================================

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_PROPERTIES_H
