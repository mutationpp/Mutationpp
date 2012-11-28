
#include "ThermalConductivityAlgorithm.h"
#include "Constants.h"
#include "Numerics.h"

using namespace Numerics;

/**
 * Computes the mixture thermal conductivity using the conjugate-gradient
 * algorithm to solve the linear transport system.
 */
class ThermalConductivityCG : public ThermalConductivityAlgorithm
{
public:

    ThermalConductivityCG(CollisionDB& collisions)
        : ThermalConductivityAlgorithm(collisions), 
          m_sys(collisions.nSpecies()-1), m_alpha(collisions.nSpecies()-1)
    { }

    double conductivity(
        const double T, const double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
    
        const RealVector& mi   = m_collisions.mass();
        const RealSymMat Astar = m_collisions.Astar(T, T, nd, p_x);
        const RealSymMat Bstar = m_collisions.Bstar(T, T, nd, p_x);
        const RealSymMat nDij  = m_collisions.nDij(T, T, nd, p_x);
        const RealVector etai  = m_collisions.etai(T, T, nd, p_x);
        
        const VectorWrapper<double> x = asVector(p_x, ns);
        
        // Compute the system matrix
        double fac = 4.0 / (15.0 * KB);
        for (int i = 1; i < ns; ++i)
            m_sys(i-1,i-1) = fac * std::max(x(i) * x(i), 1.0E-99) * mi(i) / etai(i);
        
        double miij;
        double mjij;
        for (int i = 1; i < ns; ++i) {
            for (int j = i+1; j < ns; ++j) {
                miij = mi(i) / (mi(i) + mi(j));
                mjij = mi(j) / (mi(i) + mi(j));
                fac  = std::max(x(i) * x(j), 1.0E-99) / (nDij(i,j) * 25.0 * KB);
                m_sys(i-1,j-1) = fac * miij * mjij * 
                    (16.0 * Astar(i,j) + 12.0 * Bstar(i,j) - 55.0);
                m_sys(i-1,i-1) += fac * (miij * (30.0 * miij + 16.0 * mjij * 
                    Astar(i,j)) + mjij * mjij * (25.0 - 12.0 * Bstar(i,j)));
                m_sys(j-1,j-1) += fac * (mjij * (30.0 * mjij + 16.0 * miij *
                    Astar(i,j)) + miij * miij * (25.0 - 12.0 * Bstar(i,j)));
            }
        }
        
        // Now use the Conjugate-Gradient method to solve the linear system
        DiagonalPreconditioner<double> M(m_sys);
        cg(m_sys, m_alpha, x(1,ns), M, 2);
        
        // Finally compute the dot product of alpha and X
        return dot(x(1,ns), m_alpha);
    }
    
private:

    RealSymMat m_sys;
    RealVector m_alpha;
};

// Register the algorithm
ObjectProvider<ThermalConductivityCG, ThermalConductivityAlgorithm>
    lambdaCG("CG");

