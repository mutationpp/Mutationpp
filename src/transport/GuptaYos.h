#ifndef TRANSPORT_GUPTAYOS_H
#define TRANSPORT_GUPTAYOS_H

#include "Constants.h"
#include "Numerics.h"

namespace Mutation {
    namespace Transport {

/**
 * Provides the Gupta-Yos formula as a protected member function.  This class is
 * meant to be extended by classes which compute a property using the Gupta-Yos
 * formulas (ie: viscosity and translational thermal conductivity).
 */
class GuptaYos
{
protected:
    
    template <typename T, typename E1, typename E2, typename E3>
    double guptaYos(
        const Mutation::Numerics::MatExpr<T, E1>& a, 
        const Mutation::Numerics::VecExpr<T, E2>& A, 
        const Mutation::Numerics::VecExpr<T, E3>& x)
    {
        const int ns = x.size();
    
        // Compute a_av
        double sum1 = 0.0;
        double sum2 = 0.0;
        double temp;
        
        quotient = 1.0 / A;      
        
        for (int i = 0; i < ns-1; ++i) {
            for (int j = i+1; j < ns; ++j) {
                temp = quotient(i) - quotient(j);
                temp = 2.0 * x(i) * x(j) * temp * temp;
                sum1 += temp;
                sum2 += temp * a(i,j);
            }
        }
        
        double a_av = sum2 / sum1;
        
        // Next compute the value obtained by the Gupta-Yos formula
        sum1 = 0.0;
        for (int i = 0; i < ns; ++i)
            sum1 += x(i) / (A(i) + a_av);
        
        return sum1 / (1.0 - a_av * sum1);
    }

private:
    
    Mutation::Numerics::RealVector quotient;
    
}; // class GuptaYos

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_GUPTAYOS_H


