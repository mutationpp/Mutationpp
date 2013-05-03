
#include "ThermalConductivityAlgorithm.h"
#include "Constants.h"
#include "Numerics.h"
#include "CollisionDB.h"
#include "Utilities.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

template <typename Implementation>
class ThermalConductivitySPD : public ThermalConductivityAlgorithm
{
public:
    /**
     * Constructor.
     */
    ThermalConductivitySPD(ThermalConductivityAlgorithm::ARGS arguments)
        : ThermalConductivityAlgorithm(arguments), 
          m_sys(m_collisions.nSpecies()-1), m_alpha(m_collisions.nSpecies()-1)
    { }

    double thermalConductivity(
        double Th, double Te, double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
    
        const RealVector& mi    = m_collisions.mass();
        const RealSymMat& Astar = m_collisions.Astar(Th, Te, nd, p_x);
        const RealSymMat& Bstar = m_collisions.Bstar(Th, Te, nd, p_x);
        const RealSymMat& nDij  = m_collisions.nDij(Th, Te, nd, p_x);
        const RealVector& etai  = m_collisions.etai(Th, Te, nd, p_x);
        
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
        
        // Solve the linear system (this process is determined by the subclass)
        static_cast<Implementation&>(*this).solve(x(1,ns));
        
        // Finally compute the dot product of alpha and X
        return dot(x(1,ns), m_alpha);
    }
    
protected:

    RealSymMat m_sys;
    RealVector m_alpha;
};    

/**
 * Computes the mixture thermal conductivity using the conjugate-gradient
 * algorithm to solve the linear transport system.
 */
class ThermalConductivityCG : 
    public ThermalConductivitySPD<ThermalConductivityCG>
{
    using ThermalConductivitySPD<ThermalConductivityCG>::m_alpha;
    using ThermalConductivitySPD<ThermalConductivityCG>::m_sys;

public:

    ThermalConductivityCG(CollisionDB& collisions)
        : ThermalConductivitySPD<ThermalConductivityCG>::
            ThermalConductivitySPD(collisions)
    { }

    template <typename E>
    void solve(const VecExpr<double, E>& x)
    {
        DiagonalPreconditioner<double> M(m_sys);
        cg(m_sys, m_alpha, x, M, 5);
    }
};

// Register the conjugate-gradient algorithm
Config::ObjectProvider<
    ThermalConductivityCG, ThermalConductivityAlgorithm> lambdaCG("CG");

/**
 * Computes the mixture thermal conductivity using the LDL^T decomposition
 * to solve the linear transport system.
 */
class ThermalConductivityLDLT : 
    public ThermalConductivitySPD<ThermalConductivityLDLT>
{
    using ThermalConductivitySPD<ThermalConductivityLDLT>::m_alpha;
    using ThermalConductivitySPD<ThermalConductivityLDLT>::m_sys;

public:

    ThermalConductivityLDLT(CollisionDB& collisions)
        : ThermalConductivitySPD<ThermalConductivityLDLT>::
            ThermalConductivitySPD(collisions)
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

// Register the direct solver
Config::ObjectProvider<ThermalConductivityLDLT, ThermalConductivityAlgorithm> lambdaLDLT("LDLT");


    } // namespace Transport
} // namespace Mutation

