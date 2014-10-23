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

        // Compute mass averaged collision frequency 
        sum = 0.0;
        for (int is = 1; is < ns; ++is) {
            mwis = mp_thermo->speciesMw(is);
            nu = sqrt(RU*8.0*Te/(PI*mwis))*p_X[is]*ns*Q11(is);
            sum += nu/mwis;
        } 
        sum *= 2.0*mwel;

        // Return tau
        tau = 1.0/sum;
        return tau;
      }
      
    } // namespace Transfer
} // namespace Mutation
