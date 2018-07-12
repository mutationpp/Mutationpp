/**
 * @file GSIRateManagerGamma.cpp
 *
 * @brief Class which computes the chemical production rate for each species
 *        based on tha gamma model for catalysis and ablation.
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


#include "Thermodynamics.h"
#include "Transport.h" // Cross dependence to be investigated

#include "GSIReaction.h"
#include "GSIRateLaw.h"
#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"
#include "SurfaceProperties.h"
#include "WallState.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerGamma : public GSIRateManager
{
public:
    GSIRateManagerGamma(DataGSIRateManager args)
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

    ~GSIRateManagerGamma(){}

//=============================================================================

    Eigen::VectorXd computeRate()
    {
        // Get reaction rate constant
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_react_rate_const(i_r) =
                v_reactions[i_r]->getRateLaw()->
                    forwardReactionRateCoefficient(
                        m_wall_state.getWallRhoi(), m_wall_state.getWallT());
        }

        // Constant rate times densities of species
        mv_work.setZero();
        m_reactants.incrSpecies(mv_react_rate_const, mv_work);
        m_irr_products.decrSpecies(mv_react_rate_const, mv_work);

        // Multiply by molar mass
        return mv_work.cwiseProduct(m_thermo.speciesMw().matrix());

    }

private:
    const size_t m_ns;
    const size_t m_nr;

    Eigen::VectorXd mv_react_rate_const;
    Eigen::VectorXd mv_work;

    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;

};

ObjectProvider<
    GSIRateManagerGamma, GSIRateManager>
    gsi_rate_manager_gamma("gamma");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
