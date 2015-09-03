/**
 * @file ThermalConductivityChapmannEnskog.cpp
 *
 * @brief Provides Chapmann-Enskog thermal conductivity and thermal diffusion
 * ratios using LDLT and Conjugate-Gradient linear system solvers.
 */

/*
 * Copyright 2014-2015 von Karman Institute for Fluid Dynamics (VKI)
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
#include "CollisionDB.h"
#include "Utilities.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Mutation::Utilities;
using namespace Eigen;

namespace Mutation {
    namespace Transport {

/**
 * Provides heavy particle thermal conductivity and thermal diffusion ratios
 * using the Chapmann-Enskog derivation.
 *
 * @param Solver an Eigen linear system solver
 */
template <template <typename, int> class Solver>
class ThermalConductivityChapmannEnskog : public ThermalConductivityAlgorithm
{
public:
    ThermalConductivityChapmannEnskog(CollisionDB& collisions) :
        ThermalConductivityAlgorithm(collisions),
        m_sys(MatrixXd::Zero(collisions.nHeavy(), collisions.nHeavy())),
        m_x(collisions.nHeavy()),
        m_alpha(collisions.nHeavy())
    { }

    double thermalConductivity(
        double Th, double Te, double nd, const double *const p_x)
    {
        updateAlphas(Th, Te, nd, p_x);
        return m_x.matrix().dot(m_alpha);
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
        int ik, jk;

        const ArrayXd&  mi   = m_collisions.mass();
        const ArrayXXd& nDij = m_collisions.nDij(Th, Te, nd, p_x);
        const ArrayXXd& Cij  = m_collisions.Cstij(Th, Te, nd, p_x);
        Map<const ArrayXd> xh(p_x+k, nh);

        // Compute heavy-particle thermal diffusion ratios
        m_sys.triangularView<StrictlyLower>() =
            ((xh.matrix()*xh.matrix().transpose()).array()*
            (1.2*Cij - 1.) / nDij.bottomRightCorner(nh, nh)).matrix();
        m_sys.triangularView<StrictlyUpper>() =
            m_sys.triangularView<StrictlyLower>().transpose();
        m_sys.diagonal().setZero();

        Map<ArrayXd>(p_k,ns).setZero();
        for (int j = 0; j < nh; ++j)
            Map<ArrayXd>(p_k+k,nh) += m_sys.col(j).array()*(
                m_alpha(j)*mi.tail(nh)-mi(j+k)*m_alpha.array()) /
                (mi.tail(nh) + mi(j+k));
        Map<ArrayXd>(p_k+k,nh) /= KB;

        // Do we need to compute electron contributions?
        if (k == 0)
            return;

        const ArrayXXd& Q11   = m_collisions.Q11(Th, Te, nd, p_x);
        const ArrayXXd& Q22   = m_collisions.Q22(Th, Te, nd, p_x);
        const ArrayXd&  Q12ei = m_collisions.Q12ei(Th, Te, nd, p_x);
        const ArrayXd&  Q13ei = m_collisions.Q13ei(Th, Te, nd, p_x);
        const ArrayXd&  Q14ei = m_collisions.Q14ei(Th, Te, nd, p_x);
        const ArrayXd&  Q15ei = m_collisions.Q15ei(Th, Te, nd, p_x);
        const double    Q23ee = m_collisions.Q23ee(Th, Te, nd, p_x);
        const double    Q24ee = m_collisions.Q24ee(Th, Te, nd, p_x);

        // Compute the lambdas
        static ArrayXd lam01; lam01.resize(ns);
        static ArrayXd lam02; lam02.resize(ns);

        lam01.tail(nh) = Te/Th*xh*(3.*Q12ei-2.5*Q11.col(0)).tail(nh);
        lam01(0) = -Th/Te*lam01.tail(nh).sum();
        lam02.tail(nh) = Te/Th*xh*(10.5*Q12ei-6.*Q13ei-35./8.*Q11.col(0));
        lam02(0) = -Th/Te*lam02.tail(nh).sum();

        double lam11 = p_x[0]*SQRT2*Q22(0) +
            (xh*(6.25*Q11.col(0)-15.*Q12ei+12.*Q13ei).tail(nh)).sum();
        double lam12 = p_x[0]*SQRT2*(1.75*Q22(0)-2.*Q23ee) +
            (xh*(175./16.*Q11.col(0)-315./8.*Q12ei+57.*Q13ei-30.*Q14ei).tail(nh)).sum();
        double lam22 = p_x[0]*SQRT2*(77./16.*Q22(0)-7.*Q23ee+5.*Q24ee) +
            (xh*(1225./64.*Q11.col(0)-735./8.*Q12ei+399./2.*Q13ei-210.*Q14ei+90.*Q15ei).tail(nh)).sum();

         Map<ArrayXd>(p_k,ns) += (2.5*Te/Th*p_x[0]/(lam11*lam22-lam12*lam12))*
             (lam22*lam01 - lam12*lam02);
    }

private:

    void updateAlphas(double Th, double Te, double nd, const double *const p_x)
    {
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();

        const ArrayXd&  mi    = m_collisions.mass();
        const ArrayXXd& Astar = m_collisions.Astar(Th, Te, nd, p_x);
        const ArrayXXd& Bstar = m_collisions.Bstar(Th, Te, nd, p_x);
        const ArrayXXd& nDij  = m_collisions.nDij(Th, Te, nd, p_x);
        const ArrayXd&  etai  = m_collisions.etai(Th, Te, nd, p_x);

        // Compute the system matrix
        double fac = 4.0 / (15.0 * KB);
        int k = ns - nh, ik, jk;
        m_x = Map<const ArrayXd>(p_x+k,nh).max(1.0e-16);
        m_sys.diagonal().array() = fac*m_x*m_x*mi.tail(nh)/etai.tail(nh);

        double miij;
        double mjij;
        for (int j = k; j < ns; ++j) {
            jk = j-k;
            for (int i = j+1; i < ns; ++i) {
                ik = i-k;
                miij = mi(i) / (mi(i) + mi(j));
                mjij = mi(j) / (mi(i) + mi(j));
                fac  = m_x(ik) * m_x(jk) / (nDij(i,j) * 25.0 * KB);
                m_sys(ik,jk) = fac * miij * mjij *
                    (16.0 * Astar(i,j) + 12.0 * Bstar(i,j) - 55.0);
                m_sys(ik,ik) += fac * (miij * (30.0 * miij + 16.0 * mjij *
                    Astar(i,j)) + mjij * mjij * (25.0 - 12.0 * Bstar(i,j)));
                m_sys(jk,jk) += fac * (mjij * (30.0 * mjij + 16.0 * miij *
                    Astar(i,j)) + miij * miij * (25.0 - 12.0 * Bstar(i,j)));
            }
        }

        // Solve the linear system using the type of system solution given in
        // template parameter
        solver.compute(m_sys);
        m_alpha = solver.solve(m_x.matrix());
    }

private:

    MatrixXd m_sys;
    ArrayXd m_x;
    VectorXd m_alpha;
    Solver<MatrixXd, Lower> solver;

}; // class ThermalConductivityChapmannEnskog

// Register the Chapmann-Enskog solution using the LDLT decomposition
Config::ObjectProvider<
    ThermalConductivityChapmannEnskog<Eigen::LDLT>, ThermalConductivityAlgorithm>
    lambda_CE_LDLT("Chapmann-Enskog_LDLT");

// Register the Chapmann-Enskog solution using the Conjugate-Gradient algorithm.
// Note the use of a proxy class CG to reduce the template arguments to 2
template <typename MatrixType, int UpLo>
class CG : public ConjugateGradient<MatrixType, UpLo> { };
Config::ObjectProvider<
    ThermalConductivityChapmannEnskog<CG>, ThermalConductivityAlgorithm>
    lambda_CE_CG("Chapmann-Enskog_CG");

    } // namespace Transport
} // namespace Mutation

