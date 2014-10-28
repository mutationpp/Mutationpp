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
         * \f[ \Omega^{CE}_{mi} = \sum c_1 e^E_{mi} \dot{\omega}_i\f]
         * 
         * with \f$ c_1 \f$ equal to 1 for non-preferential models.
         * 
         */
        double const OmegaCE::compute_source_Candler()
        {


        // Getting Vibrational Energy
        double* const p_helec = new double[p_thermo->nSpecies()];
        p_thermo->speciesHOverRT(NULL, NULL, NULL, NULL, p_helec, NULL);
        
        // Getting Production Rate
        double* const p_prod_rate = new double[p_thermo->nSpecies()];
        p_kinetics->netProductionRates(p_prod_rate);
        
        // Inner Product
        double c1 = 1.0E0;
        double sum = 0.E0; 

        for( int i_CE = 0 ; i_CE < p_thermo->nSpecies(); ++i_CE ){
            sum += p_prod_rate[i_CE]* p_helec[i_CE] *p_thermo->T()*RU/p_thermo->speciesMw(i_CE);
        }

        delete [] p_helec;
        delete [] p_prod_rate;
	
        return(c1*sum);
       
        }
      
    } // namespace Kinetics
} // namespace Mutation
