/**
 * @file ThermalConductivityCG.cpp
 *
 * @brief Provides thermal conductivity algorithms using CG and LDLT methods.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


#include "ThermalConductivityAlgorithm.h"
#include "Constants.h"
#include "Numerics.h"
#include "CollisionDB.h"
#include "Utilities.h"

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

template <typename SystemSolution>
class ThermalConductivitySPD : public ThermalConductivityAlgorithm
{
public:
    /**
     * Constructor.
     */
    ThermalConductivitySPD(ThermalConductivityAlgorithm::ARGS arguments)
        : ThermalConductivityAlgorithm(arguments), 
          m_sys(m_collisions.nHeavy()), m_alpha(m_collisions.nHeavy())
    { }

    double thermalConductivity(
        double Th, double Te, double nd, const double *const p_x)
    {
        // First solve the linear system for the alphas
        updateAlphas(Th, Te, nd, p_x);
        
        // Finally compute the dot product of alpha and X
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();
        const int k  = ns - nh;
        
        double lambda = 0.0;
        for (int i = k; i < ns; ++i)
            lambda += p_x[i]*m_alpha(i-k);
        
        return lambda;
    }
    
    void thermalDiffusionRatios(
        double Th, double Te, double nd, const double* const p_x,
        double* const p_k)
    {
        // First solve the linear system for the alphas
        updateAlphas(Th, Te, nd, p_x);
        
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();
        const int k  = ns - nh;
        
        const RealVector& mi   = m_collisions.mass();
        const RealSymMat& nDij = m_collisions.nDij(Th, Te, nd, p_x);
        const RealSymMat& Cij  = m_collisions.Cstij(Th, Te, nd, p_x);
        
        // Compute heavy-particle thermal diffusion ratios
        double fac;
        for (int i = k; i < ns; ++i) {
            p_k[i] = 0.0;
            for (int j = k; j < ns; ++j) {
                if (j == i) continue;
                fac = 0.1*p_x[i]*p_x[j]*(12.0*Cij(i-k,j-k)-10.0)/
                    (KB*nDij(i,j)*(mi(i)+mi(j)));
                p_k[i] += fac*m_alpha(j-k)*mi(i);
                p_k[i] -= fac*m_alpha(i-k)*mi(j);
            }
        }
        
        // Do we need to compute electron contributions?
        if (k == 0)
            return;
        
        p_k[0] = 0.0;
        if (p_x[0] < 1.0e-16)
            return;
        
        // Next compute the Lambda terms for the electron thermal diffusion
        // ratios
        const RealSymMat& Q11   = m_collisions.Q11(Th, Te, nd, p_x);
        const RealSymMat& Q22   = m_collisions.Q22(Th, Te, nd, p_x);
        const RealSymMat& B     = m_collisions.Bstar(Th, Te, nd, p_x);
        const RealVector& Q12ei = m_collisions.Q12ei(Th, Te, nd, p_x);
        const RealVector& Q13ei = m_collisions.Q13ei(Th, Te, nd, p_x);
        const RealVector& Q14ei = m_collisions.Q14ei(Th, Te, nd, p_x);
        const RealVector& Q15ei = m_collisions.Q15ei(Th, Te, nd, p_x);
        const double      Q23ee = m_collisions.Q23ee(Th, Te, nd, p_x);
        const double      Q24ee = m_collisions.Q24ee(Th, Te, nd, p_x);
        const double      me    = mi(0);
        
        // Compute the lambdas
        static RealVector lam01(ns); lam01 = 0.0;
        static RealVector lam02(ns); lam02 = 0.0;
        
        double lam11ee = 0.0;
        double lam12ee = 0.0;
        double lam22ee = 0.0;
        
        for (int i = 1; i < ns; ++i) {
            lam01(i) = -p_x[i]*(2.5*Q11(i)-3.0*Q12ei(i));
            lam01(0) -= lam01(i);
            lam02(i) = -p_x[i]*(35.0/8.0*Q11(i)-10.5*Q12ei(i)+6.0*Q13ei(i));
            lam02(0) -= lam02(i);
            lam11ee += p_x[i]*Q11(i)*(25.0/4.0-3.0*B(i));
            lam12ee += p_x[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                Q14ei(i));
            lam22ee += p_x[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                210.0*Q14ei(i)+90.0*Q15ei(i));
        }
        
        
        fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*p_x[0];
        lam01   *= fac;
        lam02   *= fac;
        lam11ee = fac*(lam11ee + SQRT2*p_x[0]*Q22(0));
        lam12ee = fac*(lam12ee + SQRT2*p_x[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
        lam22ee = fac*(lam22ee + SQRT2*p_x[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
        
        fac = 1.0/(lam11ee*lam22ee-lam12ee*lam12ee);
        for (int i = 0; i < ns; ++i)
            p_k[i] += 2.5*p_x[0]*(lam01(i)*lam22ee - lam02(i)*lam12ee)*fac;
    }
    
protected:

    RealSymMat m_sys;
    RealVector m_alpha;
    
private:

    void updateAlphas(double Th, double Te, double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();
    
        const RealVector& mi    = m_collisions.mass();
        const RealSymMat& Astar = m_collisions.Astar(Th, Te, nd, p_x);
        const RealSymMat& Bstar = m_collisions.Bstar(Th, Te, nd, p_x);
        const RealSymMat& nDij  = m_collisions.nDij(Th, Te, nd, p_x);
        const RealVector& etai  = m_collisions.etai(Th, Te, nd, p_x);
        
        const VectorWrapper<double> x = asVector(p_x, ns);
        
        // Compute the system matrix
        double fac = 4.0 / (15.0 * KB);
        int k = ns - nh;
        for (int i = k; i < ns; ++i)
            m_sys(i-k,i-k) = fac * std::max(x(i) * x(i), 1.0E-99) * mi(i) / etai(i);
        
        double miij;
        double mjij;
        for (int i = k; i < ns; ++i) {
            for (int j = i+1; j < ns; ++j) {
                miij = mi(i) / (mi(i) + mi(j));
                mjij = mi(j) / (mi(i) + mi(j));
                fac  = std::max(x(i) * x(j), 1.0E-99) / (nDij(i,j) * 25.0 * KB);
                m_sys(i-k,j-k) = fac * miij * mjij * 
                    (16.0 * Astar(i,j) + 12.0 * Bstar(i,j) - 55.0);
                m_sys(i-k,i-k) += fac * (miij * (30.0 * miij + 16.0 * mjij * 
                    Astar(i,j)) + mjij * mjij * (25.0 - 12.0 * Bstar(i,j)));
                m_sys(j-k,j-k) += fac * (mjij * (30.0 * mjij + 16.0 * miij *
                    Astar(i,j)) + miij * miij * (25.0 - 12.0 * Bstar(i,j)));
            }
        }
        
        // Solve the linear system using the type of system solution given in
        // template parameter
        static SystemSolution sol;
        sol(m_sys, x(k,ns), m_alpha);
    }
};

/**
 * Implementation of the conjugate gradient solution algorithm.
 */
class CGSolution
{
public:
    template <typename M, typename V>
    void operator () (
        const SymMatExpr<double, M>& A, const VecExpr<double, V>& b,
        Vector<double>& x)
    {
        DiagonalPreconditioner<double> P(A);
        cg(A, x, b, P, 5);
    }
};

// Register the conjugate-gradient algorithm
Config::ObjectProvider<
    ThermalConductivitySPD<CGSolution>, ThermalConductivityAlgorithm> cg("CG");

/**
 * Implementation of the LDL^T solution algorithm.
 */
class LDLTSolution
{
public:
    template <typename M, typename V>
    void operator () (
        const SymMatExpr<double, M>& A, const VecExpr<double, V>& b,
        Vector<double>& x)
    {
        static LDLT<double> ldlt;
        ldlt.setMatrix(A);
        ldlt.solve(x, b);
    }
};

// Register the LDL^T algorithm
Config::ObjectProvider<
    ThermalConductivitySPD<LDLTSolution>, ThermalConductivityAlgorithm> ldlt("LDLT");

    } // namespace Transport
} // namespace Mutation

