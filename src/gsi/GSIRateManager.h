/**
 * @file GSIRateManager.h
 *
 * @brief Declaration of abstract GSIRateManager class.
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


#ifndef GSI_RATE_MANAGER_H
#define GSI_RATE_MANAGER_H

#include <string>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction;
class SurfaceState;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the GSIRateManager class.
 */
struct DataGSIRateManager {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const SurfaceState& s_surf_state;
    const std::vector<GSIReaction*>& s_reactions;
};

//==============================================================================

/*
 * Abstract class responsible for computing the chemical production rate for
 * Gas-Surface Interaction reactions. Its implementation is different for
 * macroscopic and microscopic approaches.
 */
class GSIRateManager
{
public:
	/**
	 * Required for self registering different rate laws.
	 */
    typedef DataGSIRateManager ARGS;

	/// Returns name of this type.
	static std::string typeName() { return "GSIRateManager"; }

//==============================================================================

    /**
     * Constructor of the abstract base class, taking as input a structure
     * of type DataGSIRateManager.
     */
    GSIRateManager(ARGS args)
       : m_thermo(args.s_thermo ),
         m_surf_state(args.s_surf_state ),
         v_reactions(args.s_reactions ) { }

//==============================================================================

    /**
     * Destructor
     */
    virtual ~GSIRateManager(){ }

//==============================================================================

    /**
     * This purely virtual function returns the surface reaction rates
     * according to the selected gas-surface interaction model.
     *
     * @return Vector of the gas phase species production rates in kg/m^2-s
     */
    virtual Eigen::VectorXd computeRates() = 0;

//==============================================================================
    /**
     * This purely virtual function returns the surface production rates
     * according to the selected gas surface interaction model for each
     * reaction.
     *
     * @return Vector of size number of reactions containing
     * the reaction rates in kg/m^2-s
     */
    virtual Eigen::VectorXd computeRatesPerReaction() = 0;

//==============================================================================
    /**
     * Purely virtual function which return the number of surface reactions.
     */
    virtual int nSurfaceReactions() = 0;

//==============================================================================
protected:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceState& m_surf_state;
    const std::vector<GSIReaction*>& v_reactions;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSI_RATE_MANAGER_H
