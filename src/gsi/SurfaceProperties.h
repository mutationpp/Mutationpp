/**
 * @file SurfaceProperties.h
 *
 * @brief Declaration of SurfaceProperties class.
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


#ifndef SURFACE_PROPERTIES_H
#define SURFACE_PROPERTIES_H

#include <eigen3/Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Structure which stores the necessary inputs for the SurfaceProperties class.
 */
struct DataSurfaceProperties
{
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Utilities::IO::XmlElement& s_node_surf_props;
};

//==============================================================================

class SurfaceProperties
{
public:
	/**
	 *
	 */
    typedef const DataSurfaceProperties& ARGS;

	/// Returns name of this type.
	static std::string typeName() { return "SurfaceProperties"; }
//==============================================================================

    /**
     *
     */
    SurfaceProperties(ARGS args){ }

//==============================================================================

    /**
     *
     */
    virtual ~SurfaceProperties(){ }

//==============================================================================

    virtual int globalSpeciesIndex(const std::string& str_sp) const {
        return -1;
    };

//==============================================================================

    virtual size_t nWallSpecies() const { return 0; }

//==============================================================================

    virtual size_t nSiteCategories() const { return 0; }

//==============================================================================

    virtual double nTotalSites() const { return 0.; }

//==============================================================================
    virtual int globalIdxtoGasIdx(const size_t& id_sp_global) const {
        return -1;
    }

//==============================================================================

    virtual int globalSurfaceCompositionSpeciesIndex(
        const std::string& str_sp) const { return -1; };

//==============================================================================

    virtual double fracSite(const int& i_site) const { return 0.; }

//==============================================================================

    virtual int nSpeciesinSite(const int& i_site) const { return 0; };

//==============================================================================

    // Functions important for surfaces in equilibrium and Bprime
    virtual bool surfaceConstraint() const {return false;}

//==============================================================================

    virtual std::string surfaceSpecies() const {return "";}

//==============================================================================

    virtual double BprimePyro() const {return 0;}

//==============================================================================

    virtual void getSurfaceCompositions(
        Eigen::VectorXd & lv_pyrolisysComposition,
        Eigen::VectorXd & lv_charComposition) const {}

};

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // SURFACE_PROPERTIES_H
