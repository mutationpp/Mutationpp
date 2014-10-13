#include "Thermodynamics.h"
#include "MillikanWhite.h"
#include "TransferModel.h"
#include <cmath>


namespace Mutation {
    namespace Transfer {
      
        double const OmegaCV::compute_source_Candler()
        {
      
        /**
         * Non-preferential Model according to Candler with
         * 
         * \f[ \Omega^{CV}_{mi} = \Sumc_1 e^V_{mi} \dot{\omega}_i\f]
         * 
         * with \f$ c_1 \f$ equal to 1 for non-preferential models. For more information:
         * 
         * G. V. Candler, MacCormack, Computation of weakly ionized hypersonic flows in thermochemical nonequilibrium, Journal of
         * Thermophysics and Heat Transfer, 1991, 5(11):266
         * 
         */

        // Getting Vibrational Energy
        double* const p_enthalpy = new double[p_thermo->nSpecies()];
        double* const p_vib_enthalpy = new double[p_thermo->nSpecies()];
	double mp_work[p_thermo->nSpecies()];
        p_thermo->speciesHOverRT(p_enthalpy, mp_work, NULL, p_vib_enthalpy, NULL, NULL);
        
        // Getting Production Rate
        double* const p_prod_rate = new double[p_thermo->nSpecies()];
        p_kinetics->netProductionRates(p_prod_rate);
        
        // Inner Product

        double c1 = 1.0E0;
        double sum = 0.E0; 

        for( int i_CV = 0 ; i_CV < p_thermo->nSpecies(); ++i_CV ){
            sum += p_prod_rate[i_CV]* p_vib_enthalpy[i_CV] *p_thermo->T()*RU/p_thermo->speciesMw(i_CV);
        }

        delete [] p_enthalpy;
        delete [] p_vib_enthalpy;
        delete [] p_prod_rate;
	
        return(c1*sum);
       
        }
      
    } // namespace Kinetics
} // namespace Mutation