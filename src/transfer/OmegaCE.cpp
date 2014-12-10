/**
 * @file OmegaCE.cpp
 *
 * @brief Implementation of OmegaCE.
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
      
        /**
         * Non-preferential Model with
         * 
         * \f[ \Omega^{CE} = e^T_{e} \dot{\omega}_e\f]
         * 
         * with \f$ c_1 \f$ equal to 1 for non-preferential models.
         * 
         */
        /*double const OmegaCE::compute_source_Candler()
        {

         // Getting Electron Translational Energy
         mp_thermo->speciesHOverRT(NULL, mp_wrk1, NULL, NULL, NULL, NULL);

         // Getting Production Rate
         mp_kinetics->netProductionRates(mp_wrk2);

         // Inner Product
         double c1 = 1.0E0;
         double sum = 0.E0; 

             sum += mp_wrk1[0]*mp_wrk2[0]*mp_thermo->T()*RU/mp_thermo->speciesMw(0);

         return(c1*sum);
         }*/

         double const OmegaCE::compute_source_Candler()
         {
             mp_kinetics->netProductionRates(mp_wrk1);
             return mp_wrk1[0] * 1.5*RU*mp_thermo->Te()/mp_thermo->speciesMw(0);
         }
      
    } // namespace Kinetics
} // namespace Mutation
