/**
 * @file SurfacePropertiesNull.cpp
 *
 * @brief SurfaceProperties class when no surface properties are required
 *        to be stored.
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
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Class invoked when no surface properties are required, like in the case
 * of a purely catalytic surface modeled with the gamma model.
 */
class SurfacePropertiesNull : public SurfaceProperties
{
public:
/**
 * Default constructor
 */
	SurfacePropertiesNull(ARGS args) : SurfaceProperties(args) { }

//==============================================================================
/*
 * Default destructor
 */
    ~SurfacePropertiesNull(){ }
};

ObjectProvider<
    SurfacePropertiesNull, SurfaceProperties>
    surface_properties_null("none");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
