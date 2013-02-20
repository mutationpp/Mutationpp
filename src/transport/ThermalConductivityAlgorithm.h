#ifndef TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
#define TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H

namespace Mutation {
    namespace Transport {

class CollisionDB;

/**
 * Abstract base class for all thermal conductivity algorithms which allows for 
 * self registration of concrete types.
 */
class ThermalConductivityAlgorithm
{
public:

    /// All ThermalConductivityAlgorithms must provide a constructor taking
    /// these arguements.
    typedef CollisionDB& ARGS;

    /**
     * Constructor.
     */
    ThermalConductivityAlgorithm(ARGS collisions)
        : m_collisions(collisions)
    { }
    
    /**
     * Destructor.
     */
    virtual ~ThermalConductivityAlgorithm() {}
    
    /**
     * Returns the mixture thermal conductivity in W/m-K using the algorithm 
     * defined in the derived type.
     */
    virtual double thermalConductivity(
        double Th, double Te, double nd, const double* const p_x) = 0;
    

protected:
    
    CollisionDB& m_collisions;

}; // class ThermalConductivityAlgorithm

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
