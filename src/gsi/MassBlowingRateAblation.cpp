/**
 * @file MassBlowingRateAblation.cpp
 *
 * @brief Class which computes the mass blowing flux produced due to
 *        ablation reactions.
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


#include <Eigen/Dense>

#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceChemistry.h"
#include "MassBlowingRate.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRateAblation : public MassBlowingRate
{
public:
//==============================================================================
    /**
     * Constructor
     */
    MassBlowingRateAblation(ARGS args)
        : mv_wrk(args.s_thermo.nSpecies()),
          m_surf_chem(args.s_surf_chem)
    {}

//==============================================================================
    /**
     * Destructor
     */
    ~MassBlowingRateAblation(){}

//==============================================================================
    /**
     * This function returns the mass blowing flux in kg/m^2-s as the sum of
     * the heterogeneous reactions production rates on the surface.
     */
    double computeBlowingFlux(){
        m_surf_chem.surfaceReactionRates(mv_wrk);
        return mv_wrk.sum();
    }

//==============================================================================
    /**
     * This function returns the mass blowing flux in kg/m^2-s as the sum of
     * the heterogeneous reactions production rates on the surface
     * given the chemical production rates for all species.
     */
    double computeBlowingFlux(const Eigen::VectorXd& v_chem_rates){
        return v_chem_rates.sum();
    }

private:
    const SurfaceChemistry& m_surf_chem;

    VectorXd mv_wrk;
};

ObjectProvider<
    MassBlowingRateAblation, MassBlowingRate>
    mass_blowing_rate_ablation("isOn");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
