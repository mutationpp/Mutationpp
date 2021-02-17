/**
 * @file GSIRateManagerPhenomenological.cpp
 *
 * @brief Class which computes the chemical production rate for each species
 *        based on phenomenological models for catalysis and ablation.
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

#include "GSIReaction.h"
#include "GSIRateLaw.h"
#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"
#include "SurfaceProperties.h"
#include "SurfaceState.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerPhenomenological : public GSIRateManager
{
public:
    GSIRateManagerPhenomenological(DataGSIRateManager args)
        : GSIRateManager(args),
		  m_ns(args.s_thermo.nSpecies()),
		  m_nr(args.s_reactions.size()),
          mv_react_rate_const(m_nr),
		  mv_work(m_ns)
    {
        for (int i_reac = 0; i_reac < m_nr; ++i_reac) {
            m_reactants.addReaction(
                i_reac, args.s_reactions[i_reac]->getReactants());
            m_irr_products.addReaction(
                i_reac, args.s_reactions[i_reac]->getProducts());
        }
    }

//=============================================================================

    ~GSIRateManagerPhenomenological(){}

//=============================================================================

    Eigen::VectorXd computeRates()
    {
        // Get reaction rate constant
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_react_rate_const(i_r) =
                v_reactions[i_r]->getRateLaw()->
                    forwardReactionRateCoefficient(
                        m_surf_state.getSurfaceRhoi(),
                        m_surf_state.getSurfaceT());
        }

        // Constant rate times densities of species
        mv_work.setZero();
        m_reactants.incrSpecies(mv_react_rate_const, mv_work);
        m_irr_products.decrSpecies(mv_react_rate_const, mv_work);

        // Multiply by molar mass
        return mv_work.cwiseProduct(m_thermo.speciesMw().matrix());
    }

//=============================================================================

    Eigen::VectorXd computeRatesPerReaction()
    {
    	// Getting the kfs with the initial conditions
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_react_rate_const(i_r) =
                v_reactions[i_r]->getRateLaw()->forwardReactionRateCoefficient(
            		m_surf_state.getSurfaceRhoi(), m_surf_state.getSurfaceT());
        }
        m_reactants.multReactions(
            m_surf_state.getSurfaceRhoi(), mv_react_rate_const);

        return mv_react_rate_const;
    }


//=============================================================================

    int nSurfaceReactions(){ return m_nr; }

//=============================================================================
private:
    const size_t m_ns;
    const size_t m_nr;

    Eigen::VectorXd mv_react_rate_const;
    Eigen::VectorXd mv_work;

    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;
};

ObjectProvider<GSIRateManagerPhenomenological, GSIRateManager>
    gsi_rate_manager_phenomenological_mass("phenomenological_mass");

ObjectProvider<GSIRateManagerPhenomenological, GSIRateManager>
    gsi_rate_manager_phenomenological_mass_energy(
        "phenomenological_mass_energy");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
