/**
 * @file MassBlowingRate.h
 *
 * @brief Declaration of MassBlowingRate class.
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


#ifndef MASS_BLOWING_RATE_H
#define MASS_BLOWING_RATE_H

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceChemistry;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the MassBlowingRate class.
 */
struct DataMassBlowingRate {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const SurfaceChemistry& s_surf_chem;
};

//==============================================================================

/**
 * Abstract class which returns the mass blowing rate for a heterogeneous
 * reactions. It is the base of a Null object (not used yet) which returns
 * zero blow flux. This is important for the case of purely catalytic or
 * oxidation reactions. In the case where ablative reactions take place,
 * a solid class ablation is available.
 */

class MassBlowingRate
{
public:
	/**
	 * Required for self registering different type of mass blowing rate.
	 */
    typedef const DataMassBlowingRate& ARGS;

	/// Returns name of this type.
	static std::string typeName() { return "MassBlowingRate"; }

    /**
     * Destructor
     */
    virtual ~MassBlowingRate(){ }

    /**
     * Purely virtual function which return the blowing flux in kg/m^2-s.
     */
    virtual double computeBlowingFlux() = 0;

    /**
     * Purely virtual function which return the blowing flux in kg/m^2-s
     * given the chemical production rates for all species.
     */
    virtual double computeBlowingFlux(
        const Eigen::VectorXd& v_chem_rates) = 0;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // MASS_BLOWING_RATE_H
