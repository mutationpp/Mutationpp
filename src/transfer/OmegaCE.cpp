#include "Thermodynamics.h"
#include "MillikanWhite.h"
#include "TransferModel.h"
#include <cmath>


namespace Mutation {
    namespace Transfer {
      
        double const OmegaCE::compute_source_Candler()
        {
      
        /**
         * Non-preferential Model with
         * 
         * \f[ \Omega^{CE}_{mi} = \Sumc_1 e^E_{mi} \dot{\omega}_i\f]
         * 
         * with \f$ c_1 \f$ equal to 1 for non-preferential models.
         * 
         */

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
