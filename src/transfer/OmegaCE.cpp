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
         * \f[ \Omega^{CE} = e^T_{e} \dot{\omega}_e\f]
         * 
         * with \f$ c_1 \f$ equal to 1 for non-preferential models.
         * 
         */

         // Getting Electron Translational Energy
         mp_thermo->speciesHOverRT(NULL, mp_wrk1, NULL, NULL, NULL, NULL);

         // Getting Production Rate
         mp_kinetics->netProductionRates(mp_wrk2);

         // Inner Product
         double c1 = 1.0E0;
         double sum = 0.E0; 

             sum += mp_wrk1[0]*mp_wrk2[0]*mp_thermo->T()*RU/mp_thermo->speciesMw(0);

         return(c1*sum);
         }
      
    } // namespace Kinetics
} // namespace Mutation
