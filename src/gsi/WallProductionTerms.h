/**
 * @file WallProductionTerms.h
 *
 * @brief Declaration of abstract WallProductionTerms class.
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


#ifndef WALL_PRODUCTION_TERMS_H
#define WALL_PRODUCTION_TERMS_H

#include <eigen3/Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction;
class GSIRateManager;
class WallState;
class SurfaceProperties;
class WallProductionTerms;

/**
 * Structure which stores the necessary inputs for the
 * WallProductionTerms class.
 */
struct DataWallProductionTerms {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism;
    const Mutation::Utilities::IO::XmlElement& s_node_prod_terms;
    const SurfaceProperties& s_surf_props;
    const WallState& s_wall_state;
    std::vector<WallProductionTerms*>* sp_surf_prod;
    const double* const sp_pres;
};

//==============================================================================
/*
 *  Abstract class
 */
class WallProductionTerms
{
public:
//==============================================================================
    typedef const DataWallProductionTerms& ARGS;

//==============================================================================
    // Constructor
    WallProductionTerms(ARGS args){}

//==============================================================================
	/// Returns name of this type.
	static std::string typeName() { return "WallProductionTerms"; }

//==============================================================================
    virtual ~WallProductionTerms(){}

//==============================================================================
    virtual void productionRate(Eigen::VectorXd& lv_mass_prod_rate) = 0;

//==============================================================================
    virtual const std::string& getWallProductionTermTag() const = 0;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // WALL_PRODUCTION_TERMS_H
