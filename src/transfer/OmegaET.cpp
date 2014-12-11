/**
 * @file OmegaET.cpp
 *
 * @brief Implementation of OmegaET.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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
#include "MillikanWhite.h"
#include "TransferModel.h"
#include <cmath>


namespace Mutation {
    namespace Transfer {
      
      double const OmegaET::compute_tau_ET(){

        const double T = mp_thermo->T();
        const double Te = mp_thermo->Te();
        const double nd = mp_thermo->numberDensity();    
        const double ns = mp_thermo->nSpecies();
        const double *const p_X = mp_thermo->X();
        const double mwel = mp_thermo->speciesMw(0);
        double mwis, nu, sum, tau; 

        // Collisional integrals
        const Numerics::RealSymMat& Q11 = mp_collisions->Q11(T, Te, nd, p_X);

        // Electron velocity
        double ve = sqrt(RU*8.0*Te/(PI*mwel));

        // Compute mass averaged collision frequency 
        sum = 0.0;
        for (int is = 1; is < ns; ++is) {
            mwis = mp_thermo->speciesMw(is);
            nu = ve * p_X[is] * nd * Q11(is);
            sum += nu/mwis;
        }
        sum *= 2.0*mwel;

        // Return tau
        tau = 1.0/sum;
        return tau;
      }
      
    } // namespace Transfer
} // namespace Mutation
