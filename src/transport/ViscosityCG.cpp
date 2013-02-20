/**
 * @file ViscosityCG.cpp
 */

#include "ViscosityAlgorithm.h"
#include "Numerics.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

template <typename Implementation>
class ViscositySPD : public ViscosityAlgorithm
{
public:

    ViscositySPD(CollisionDB& collisions)
        : ViscosityAlgorithm(collisions), m_sys(collisions.nSpecies()-1),
          m_alpha(collisions.nSpecies()-1)
    { }
    
    double viscosity(
        const double T, const double nd, const double *const p_x) 
    {
        const int ns = m_collisions.nSpecies();
        
        const RealVector& mi    = m_collisions.mass();
        const RealSymMat& Astar = m_collisions.Astar(T, T, nd, p_x);
        const RealSymMat& nDij  = m_collisions.nDij(T, T, nd, p_x);
        const RealVector& etai  = m_collisions.etai(T, T, nd, p_x);
        const VectorWrapper<double> x = asVector(p_x, ns);
        
        double fac;
        
        // Form the symmetric positive definite system matrix
        for (int i = 1; i < ns; ++i)
            m_sys(i-1,i-1) = std::max(x(i) * x(i), 1.0E-99) / etai(i);
            
        for (int i = 1; i < ns; ++i) {
            for (int j = i+1; j < ns; ++j) {
                fac = std::max(x(i) * x(j), 1.0E-99) / 
                    (nDij(i,j) * (mi(i) + mi(j)));
                m_sys(i-1,j-1) =  fac * (1.2 * Astar(i,j) - 2.0);
                m_sys(i-1,i-1) += fac * (1.2 * mi(j) / mi(i) * Astar(i,j) + 2.0);
                m_sys(j-1,j-1) += fac * (1.2 * mi(i) / mi(j) * Astar(i,j) + 2.0);
            }
        }
        
        // Solve the linear system (this process is determined by the subclass)
        static_cast<Implementation&>(*this).solve(x(1,ns));
        
        // Finally compute the dot product of alpha and X
        return dot(x(1,ns), m_alpha);
    }

protected:

    RealSymMat m_sys;
    RealVector m_alpha;     
    
}; // ViscosityCG


/**
 * Implements the conjugate gradient algorithm for solving the viscosity 
 * transport system.
 */
class ViscosityCG : public ViscositySPD<ViscosityCG>
{

    using ViscositySPD<ViscosityCG>::m_alpha;
    using ViscositySPD<ViscosityCG>::m_sys;

public:

    ViscosityCG(CollisionDB& collisions)
        : ViscositySPD<ViscosityCG>::ViscositySPD(collisions)
    { }

    template <typename E>
    void solve(const VecExpr<double, E>& x)
    {
        DiagonalPreconditioner<double> M(m_sys);
        cg(m_sys, m_alpha, x, M, 1);
    }
};

// Register this algorithm with the other ViscosityAlgorithm types
Config::ObjectProvider<ViscosityCG, ViscosityAlgorithm> viscosityCG("CG");


/**
 * Implements the LDL^T solution method for solving the viscosity transport
 * system.
 */
class ViscosityLDLT : public ViscositySPD<ViscosityLDLT>
{
    using ViscositySPD<ViscosityLDLT>::m_alpha;
    using ViscositySPD<ViscosityLDLT>::m_sys;

public:

    ViscosityLDLT(CollisionDB& collisions)
        : ViscositySPD<ViscosityLDLT>::ViscositySPD(collisions)
    { }

    template <typename E>
    void solve(const VecExpr<double, E>& x)
    {
        m_ldlt.setMatrix(m_sys);
        m_ldlt.solve(m_alpha, x);
    }

private:

    LDLT<double> m_ldlt;
};

// Register this algorithm with the other ViscosityAlgorithm types
Config::ObjectProvider<ViscosityLDLT, ViscosityAlgorithm> viscosityLDLT("LDLT");


    } // namespace Transport
} // namespace Mutation

