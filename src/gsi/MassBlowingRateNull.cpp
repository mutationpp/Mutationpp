/**
 * @file MassBlowingRateNull.cpp
 *
 * @brief Class which set the mass blowing flux equal to zero when
 *        no ablation reactions are present.
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

#include "MassBlowingRate.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateNull : public MassBlowingRate
{
public:
//==============================================================================
	/**
	 * Constructor.
	 */
	MassBlowingRateNull(ARGS args){ }

//==============================================================================
	/**
	 * Destructor.
	 */
	~MassBlowingRateNull(){ }

//==============================================================================
	/**
	 * Returns mass blowing flux equal to zero.
	 */
	double computeBlowingFlux(){ return 0.0; }

//==============================================================================
	/**
	 * Returns mass blowing flux equal to zero.
	 */
    double computeBlowingFlux(const Eigen::VectorXd& v_chem_rates){
		return 0.;
	}
};

ObjectProvider<
    MassBlowingRateNull, MassBlowingRate>
    mass_blowing_rate_null("zero");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
