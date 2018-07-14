/**
 * @file MassBlowingRateAblation.cpp
 *
 * @brief Class which computes the mass blowing flux produced due to
 *        ablation reactions.
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


#include <eigen3/Eigen/Dense>

#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "MassBlowingRate.h"
#include "WallProductionTerms.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateAblation : public MassBlowingRate
{
public:
    MassBlowingRateAblation(ARGS args)
        : m_ns(args.s_thermo.nSpecies()),
          mv_wall_prod_rates(m_ns),
          mv_wrk(m_ns+args.s_thermo.nEnergyEqns())
    {
    	std::string prod_term_tag;
        for (size_t i_term = 0;
             i_term < args.vs_wall_productions_terms.size();
             ++i_term) {
            prod_term_tag = args.vs_wall_productions_terms[i_term]->
                                getWallProductionTermTag();

            if (prod_term_tag.compare("surface_chemistry") == 0) {
                mv_wall_prod_terms.push_back(
                    args.vs_wall_productions_terms[i_term]);
            }
            if(prod_term_tag.compare("pyrolysis") == 0) {
                mv_wall_prod_terms.push_back(
                    args.vs_wall_productions_terms[i_term]);
            }
        }

    }

//==============================================================================
    /**
     * Destructor
     */
    ~MassBlowingRateAblation(){}

//==============================================================================
    /**
     * This function returns the mass blowing flux in kg/m^2-s as the sum of
     * the heterogeneous reactions.
     */
    double computeBlowingFlux(){

        mv_wall_prod_rates.setZero();
        mv_wrk.setZero();

        for(size_t i_term = 0; i_term < mv_wall_prod_terms.size(); ++i_term)
        {
        	mv_wall_prod_terms[i_term]->productionRate(mv_wrk);
            mv_wall_prod_rates += mv_wrk.head(m_ns);
        }

        return mv_wall_prod_rates.sum();
    }

private:
    const size_t m_ns;

    std::vector<WallProductionTerms*> mv_wall_prod_terms;
    VectorXd mv_wall_prod_rates;
    VectorXd mv_wrk;
};

ObjectProvider<
    MassBlowingRateAblation, MassBlowingRate>
    mass_blowing_rate_ablation("isOn");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
