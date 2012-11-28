
#include "ViscosityAlgorithm.h"
#include "GuptaYos.h"
#include "Numerics.h"

using namespace Numerics;

/**
 * Uses the Gupta-Yos formula to compute viscosity.
 */
class ViscosityGuptaYos : public ViscosityAlgorithm, public GuptaYos
{
public:

    ViscosityGuptaYos(CollisionDB& collisions) 
        : ViscosityAlgorithm(collisions), a(collisions.nSpecies()-1), 
          B(collisions.nSpecies()-1), A(collisions.nSpecies()-1)
    { }

    /**
     * Returns the viscosity of the mixture in Pa-s.
     */
    double viscosity(const double T, const double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
        
        const RealSymMat& nDij  = m_collisions.nDij(T, T, nd, p_x);
        const RealSymMat& Astar = m_collisions.Astar(T, T, nd, p_x);
        const RealVector& mass  = m_collisions.mass();
        
        // Compute the symmetric matrix a and vector A
        for (int i = 1; i < ns; ++i) {
            for (int j = i; j < ns; ++j) {
                a(i-1,j-1) = (2.0 - 1.2 * Astar(i,j)) / 
                    ((mass(i) + mass(j)) * nDij(i,j));
                B(i-1,j-1) = Astar(i,j) / nDij(i,j);
            }
        }
                
        // Leave B as symmetric and postpone dividing by m(i) until after B*x
        // which uses fewer divides and allows for sym-matrix-vector product
        A = B * asVector(p_x, ns)(1,ns);
        A = 1.2 * A / mass(1,ns);
            
        // Now compute the viscosity using Gupta-Yos
        return guptaYos(a, A, asVector(p_x, ns)(1,ns));
    }
    
private:
    
    RealSymMat a;
    RealSymMat B;
    RealVector A;
};

ObjectProvider<ViscosityGuptaYos, ViscosityAlgorithm> 
    viscosityGuptaYos("Gupta-Yos");
    
    

