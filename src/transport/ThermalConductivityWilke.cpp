
#include "ThermalConductivityAlgorithm.h"
#include "Wilke.h"
#include "Constants.h"
#include "Numerics.h"

using namespace Numerics;

/**
 * Computes the translational thermal conductivity of a mixture using the Wilke 
 * mixture rule and Eucken's relation.
 *
 * @todo I think I don't actually account for Eucken's correction here...
 */
class ThermalConductivityWilke 
    : public ThermalConductivityAlgorithm, public Wilke
{
public:
 
    ThermalConductivityWilke(CollisionDB& collisions)
        : ThermalConductivityAlgorithm(collisions)
    { }
    
    double conductivity(
        const double T, const double nd, const double *const p_x) 
    {
        const size_t ns = m_collisions.nSpecies();
        const RealVector Mw = m_collisions.mass() * NA;
        
        return wilke(
            3.75 * RU * m_collisions.etai(T, T, nd, p_x) / Mw, 
            m_collisions.mass(), asVector(p_x, ns));
    }
};

// Register the algorithm
ObjectProvider<ThermalConductivityWilke, ThermalConductivityAlgorithm>
    lambdaWilke("Wilke");

