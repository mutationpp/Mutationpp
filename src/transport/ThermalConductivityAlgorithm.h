#ifndef TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
#define TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H

#include "AutoRegistration.h"
#include "CollisionDB.h"

/**
 * Abstract base class for all thermal conductivity algorithms which allows for 
 * self registration of concrete types.
 */
class ThermalConductivityAlgorithm
{
public:

    typedef CollisionDB& ARGS;
    typedef Utilities::Provider<ThermalConductivityAlgorithm> PROVIDER;

    ThermalConductivityAlgorithm(ARGS collisions)
        : m_collisions(collisions)
    { }
    
    virtual ~ThermalConductivityAlgorithm() { }
    
    /**
     * Returns the mixture thermal conductivity.s
     */
    virtual double conductivity(
        const double T, const double nd, const double *const p_x) = 0;

protected:
    
    CollisionDB& m_collisions;

}; // 

#endif // TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
