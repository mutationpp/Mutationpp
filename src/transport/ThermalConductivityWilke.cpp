
#include "ThermalConductivityAlgorithm.h"
#include "Wilke.h"
#include "Constants.h"
#include "Numerics.h"
#include "CollisionDB.h"
#include "AutoRegistration.h"

using namespace Numerics;

/**
 * Computes the translational thermal conductivity of a mixture using the Wilke 
 * mixture rule and Eucken's relation.
 */
class ThermalConductivityWilke 
    : public ThermalConductivityAlgorithm, public Wilke
{
public:
 
    ThermalConductivityWilke(ThermalConductivityAlgorithm::ARGS arguments)
        : ThermalConductivityAlgorithm(arguments)
    { }
    
    double elasticThermalConductivity(
        double Th, double Te, double nd, const double* const p_x) 
    {
        const size_t ns = m_collisions.nSpecies();
        const RealVector Mw = m_collisions.mass() * NA;
        
        return wilke(
            3.75 * RU * m_collisions.etai(Th, Te, nd, p_x) / Mw, 
            m_collisions.mass(), asVector(p_x, ns));
    }
};

// Register the algorithm
Utilities::ObjectProvider<ThermalConductivityWilke, ThermalConductivityAlgorithm>
    lambdaWilke("Wilke");

