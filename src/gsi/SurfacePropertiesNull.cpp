/**
 * @file SurfacePropertiesNull.cpp
 *
 * @brief SurfaceProperties class when no surface properties are required
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


#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfacePropertiesNull : public SurfaceProperties
{
public:
	SurfacePropertiesNull(ARGS args) : SurfaceProperties(args) { }

	//==============================================================================

    ~SurfacePropertiesNull(){ }

	//==============================================================================

    int speciesIndexWall(const std::string& str_sp) const { return -1; }

    //==============================================================================

    int nSpeciesWall() const { return 0; }

	//==============================================================================

    int nSites() const { return 0; }

    //==============================================================================

    double nTotalSites() const { return 0.0; }

	//==============================================================================

    double fracSite(const int& i_site) const { return 0.0; }

    //==============================================================================

    int nSpeciesSite(const int& i_site) const { return 0; }
};

ObjectProvider<
    SurfacePropertiesNull, SurfaceProperties>
    surface_properties_gamma_null("gamma");

    } // namespace GasSurfaceInteraction
} // namespace Mutation 
