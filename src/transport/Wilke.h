#ifndef TRANSPORT_WILKE_H
#define TRANSPORT_WILKE_H

#include "Numerics.h"

/**
 * Base class for class that compute transport properties using the Wilke 
 * formulas.
 */
class Wilke
{
public:
    
    Wilke() { }

protected:

    template <typename E1, typename E2, typename E3>
    double wilke(
        const Numerics::VecExpr<double, E1>& vals, 
        const Numerics::VecExpr<double, E2>& mass,
        const Numerics::VecExpr<double, E3>& x)
    {
        const int ns = vals.size();

        double average = 0.0;
        double sum;
        double ratio;
        double temp;
        
        for (int i = 1; i < ns; ++i) {
            sum = 0.0;
            
            for (int j = 1; j < ns; ++j) {
                if (i == j)
                    sum += x(j);
                else {
                    ratio = mass(i) / mass(j);
                    temp = 1.0 + std::sqrt(vals(i) / vals(j) / std::sqrt(ratio));
                    sum += x(j) * temp * temp / std::sqrt(8.0 * (1.0 + ratio));
                }
            }
            
            average += x(i) * vals(i) / sum;            
        }
        
        return average;
    }
    
}; // class Wilke


#endif // TRANSPORT_WILKE_H
