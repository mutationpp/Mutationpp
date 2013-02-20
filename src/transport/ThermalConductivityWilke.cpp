
#include "ThermalConductivityAlgorithm.h"
#include "Wilke.h"
#include "Constants.h"
#include "Numerics.h"
#include "CollisionDB.h"
#include "Utilities.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

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
    
    double thermalConductivity(
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
Config::ObjectProvider<ThermalConductivityWilke, ThermalConductivityAlgorithm>
    lambdaWilke("Wilke");

    } // namespace Transport
} // namespace Mutation

