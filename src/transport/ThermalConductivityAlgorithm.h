#ifndef TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
#define TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H

#include <utility>

class Thermodynamics;
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
    typedef std::pair<Thermodynamics&, CollisionDB&> ARGS;

    /**
     * Constructor.
     */
    ThermalConductivityAlgorithm(ARGS arguments);
    
    /**
     * Destructor.
     */
    virtual ~ThermalConductivityAlgorithm();
    
    /**
     * Returns the mixture thermal conductivity in W/m-K using the algorithm 
     * defined in the derived type and applying the Euken corrections for 
     * internal energy.
     */
    virtual double thermalConductivity(
        double Th, double Te, double Tr, double Tv, double Tel, double n, 
        const double* const p_x);

protected:
    
    /**
     * Computes the thermal conductivity using the algorithm implemented in the
     * derived type.  Note this should not include the application of Euken's
     * correction for inelastic collisions.
     */
    virtual double elasticThermalConductivity(
        double Th, double Te, double nd, const double* const p_x) = 0;
    

protected:
    
    Thermodynamics& m_thermo;
    CollisionDB& m_collisions;

private:
    
    double* mp_cvr;
    double* mp_cvv;
    double* mp_cvel;

}; // class ThermalConductivityAlgorithm

#endif // TRANSPORT_THERMAL_CONDUCTIVITY_ALGORITHM_H
