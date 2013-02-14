
#include "ViscosityAlgorithm.h"
#include "Wilke.h"
#include "Numerics.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

/**
 * Computes the viscosity of a mixture in Pa-s using the Wilke mixture rule.
 */
class ViscosityWilke : public ViscosityAlgorithm, public Wilke
{
public:

    ViscosityWilke(CollisionDB& collisions) 
        : ViscosityAlgorithm(collisions)
    { }

    /**
     * Returns the viscosity of the mixture in Pa-s.
     */
    double viscosity(const double T, const double nd, const double *const p_x) 
    {
        const size_t ns = m_collisions.nSpecies();
        
        return wilke(
            m_collisions.etai(T, T, nd, p_x), m_collisions.mass(), 
            asVector(p_x, ns));
    }
       
};

Config::ObjectProvider<ViscosityWilke, ViscosityAlgorithm> viscosityWilke("Wilke");


