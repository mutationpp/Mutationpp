/**
 * @file SurfaceChemistry.h
 *
 * @brief Class which takes care of the surface chemistry.
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


#ifndef SURFACE_CHEMISTRY_H
#define SURFACE_CHEMISTRY_H

#include <eigen3/Eigen/Dense>
#include <string>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction;
class GSIRateManager;
class SurfaceState;

//==============================================================================
/**
 *  Class which ...
 */
class SurfaceChemistry
{
public:
    // Constructor
    SurfaceChemistry(
        Mutation::Thermodynamics::Thermodynamics& thermo,
        const Mutation::Transport::Transport& transport,
        const std::string& gsi_mechanism,
        const Mutation::Utilities::IO::XmlElement& xml_surf_chem,
        const SurfaceState& surf_state);

//==============================================================================
    /**
     * Default destructor
     */
    ~SurfaceChemistry();

//==============================================================================
    /**
     * Chemical source term per kg
     */
    void surfaceReactionRates(Eigen::VectorXd& v_mass_chem_rate) const;

//==============================================================================

    void surfaceReactionRatesPerReaction(Eigen::VectorXd& v_rate_per_reaction);

//==============================================================================
    /**
     * What do I do if these are zero!
     */
    int nSurfaceReactions();

//==============================================================================

private:
    void addReaction(GSIReaction* gsi_reaction){
        mv_reaction.push_back(gsi_reaction);
    }

private:
    const size_t m_ns;

    std::string m_reaction_type;

    std::vector<GSIReaction*> mv_reaction;
    GSIRateManager* mp_rate_manager;
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_CHEMISTRY_H
