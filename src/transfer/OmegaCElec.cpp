#include "Thermodynamics.h"
#include "MillikanWhite.h"
#include "TransferModel.h"
#include <cmath>


namespace Mutation {
    namespace Transfer {
      
        double const OmegaCElec::compute_source_Candler()
        {
        /**
         * Non-preferential Model with
         * 
         * \f[ \Omega^{CElec}_{mi} = \Sumc_1 e^{Elec}_{mi} \dot{\omega}_i\f]
         * 
         * with \f$ c_1 \f$ equal to 1 for non-preferential models.
         * 
         */

         // Getting Electronic Energy
         mp_thermo->speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);

         // Getting Production Rate
         mp_kinetics->netProductionRates(mp_wrk2);

         // Inner Product
         double c1 = 1.0E0;
         double sum = 0.E0; 

         for(int i = 0 ; i < m_ns; ++i)
             sum += mp_wrk1[i]*mp_wrk2[i] *mp_thermo->T()*RU/mp_thermo->speciesMw(i);

         return(c1*sum);
         }
      
    } // namespace Kinetics
} // namespace Mutation
