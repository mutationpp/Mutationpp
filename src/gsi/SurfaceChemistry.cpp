/**
 * @file SurfaceChemistry.cpp
 *
 * @brief Implementation of the SurfaceChemistry class.
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


#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "GSIRateManager.h"
#include "GSIReaction.h"
#include "SurfaceChemistry.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

SurfaceChemistry::SurfaceChemistry(
    Mutation::Thermodynamics::Thermodynamics& thermo,
    const Mutation::Transport::Transport& transport,
    const std::string& gsi_mechanism,
    const Mutation::Utilities::IO::XmlElement& xml_surf_chem,
    const SurfaceState& surf_state)
        :   m_ns(thermo.nSpecies()),
            mp_rate_manager(NULL)
{
    // Filling in the reaction vector
    DataGSIReaction data_gsi_reaction = {
        thermo, transport, surf_state};

    for (data_gsi_reaction.s_iter_reaction = xml_surf_chem.begin();
         data_gsi_reaction.s_iter_reaction != xml_surf_chem.end();
         ++data_gsi_reaction.s_iter_reaction)
    {
        data_gsi_reaction.s_iter_reaction->
            getAttribute("type", m_reaction_type,
                "Error in the reaction input. A type attribute should "
              "be provided with the reaction type.");

        addReaction(Factory<GSIReaction>::create(
            m_reaction_type, data_gsi_reaction));
    }


    // Assigning the mp_rate_manager pointer
    DataGSIRateManager data_rate_manager = {
        thermo, surf_state, mv_reaction };
    mp_rate_manager = Factory<GSIRateManager>::create(
        gsi_mechanism, data_rate_manager);
}

//==============================================================================

SurfaceChemistry::~SurfaceChemistry()
{
    for (std::vector<GSIReaction*>::iterator iter = mv_reaction.begin();
    	 iter != mv_reaction.end(); ++iter){
    	delete (*iter);
    }
    mv_reaction.clear();

    if (mp_rate_manager != NULL){ delete mp_rate_manager; }
}

//==============================================================================

void SurfaceChemistry::surfaceReactionRates(VectorXd& v_mass_chem_rate) const {
    v_mass_chem_rate = mp_rate_manager->computeRates();
}

//==============================================================================

void SurfaceChemistry::surfaceReactionRatesPerReaction(
    Eigen::VectorXd& v_rate_per_reaction)
{
    v_rate_per_reaction = mp_rate_manager->computeRatesPerReaction();
}

//==============================================================================

int SurfaceChemistry::nSurfaceReactions() {
    return mp_rate_manager->nSurfaceReactions();
}

    } // namespace GasSurfaceInteraction
} // namespace Mutation
