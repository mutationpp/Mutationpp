/**
 * @file ViscosityCG.cpp
 *
 * @brief Provides viscosity classes using CG and LDLT solvers.
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

#include "ViscosityAlgorithm.h"
#include "Numerics.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

template <template <typename, int> class Solver>
class ViscosityChapmannEnskog : public ViscosityAlgorithm
{
public:
	ViscosityChapmannEnskog(CollisionDB& collisions) :
		ViscosityAlgorithm(collisions),
		m_sys(collisions.nHeavy(), collisions.nHeavy()),
		m_x(collisions.nHeavy()),
		m_alpha(collisions.nHeavy())
	{ }

	double viscosity(
	    const double T, const double nd, const double *const p_x)
	{
		const int ns = m_collisions.nSpecies();
		const int nh = m_collisions.nHeavy();

		const RealVector& mi    = m_collisions.mass();
		const RealSymMat& Astar = m_collisions.Astar(T, T, nd, p_x);
		const RealSymMat& nDij  = m_collisions.nDij(T, T, nd, p_x);
		const RealVector& etai  = m_collisions.etai(T, T, nd, p_x);

		// Form the symmetric positive definite system matrix
		int k = ns-nh, ik, jk;
		double fac;

		for (int i = k; i < ns; ++i) {
			ik = i-k;
			m_x(ik) = std::max(p_x[i], 1.0e-16);
			m_sys(ik,ik) = m_x(ik) * m_x(ik) / etai(i);
		}

		for (int i = k; i < ns; ++i) {
			ik = i-k;
			for (int j = i+1; j < ns; ++j) {
				jk = j-k;
				fac = m_x(ik)*m_x(jk) / (nDij(i,j) * (mi(i) + mi(j)));
				m_sys(ik,jk) =  fac * (1.2 * Astar(i,j) - 2.0);
				m_sys(ik,ik) += fac * (1.2 * mi(j) / mi(i) * Astar(i,j) + 2.0);
				m_sys(jk,jk) += fac * (1.2 * mi(i) / mi(j) * Astar(i,j) + 2.0);
			}
		}

		// Solve the linear system
		m_solver.compute(m_sys);
		m_alpha = m_solver.solve(m_x);

		// Finally compute the dot product of alpha and X
		return m_x.dot(m_alpha);
	}

private:

	Eigen::MatrixXd m_sys;
	Eigen::VectorXd m_x;
	Eigen::VectorXd m_alpha;
	Solver<Eigen::MatrixXd, Eigen::Upper> m_solver;
};

// Register the Chapmann-Enskog solution using the LDLT decomposition
Config::ObjectProvider<
	ViscosityChapmannEnskog<Eigen::LDLT>, ViscosityAlgorithm>
	visc_CE_LDLT("Chapmann-Enskog_LDLT");

// Register the Chapmann-Enskog solution using the Conjugate-Gradient algorithm.
// Note the use of a proxy class CG to reduce the template arguments to 2
template <typename MatrixType, int UpLo>
class CG : public Eigen::ConjugateGradient<MatrixType, UpLo> { };
Config::ObjectProvider<
	ViscosityChapmannEnskog<CG>, ViscosityAlgorithm>
	visc_CE_CG("Chapmann-Enskog_CG");



//template <typename Implementation>
//class ViscositySPD : public ViscosityAlgorithm
//{
//public:
//
//    ViscositySPD(CollisionDB& collisions)
//        : ViscosityAlgorithm(collisions), m_sys(collisions.nHeavy()),
//          m_alpha(collisions.nHeavy())
//    { }
//
//    double viscosity(
//        const double T, const double nd, const double *const p_x)
//    {
//        const int ns = m_collisions.nSpecies();
//        const int nh = m_collisions.nHeavy();
//
//        const RealVector& mi    = m_collisions.mass();
//        const RealSymMat& Astar = m_collisions.Astar(T, T, nd, p_x);
//        const RealSymMat& nDij  = m_collisions.nDij(T, T, nd, p_x);
//        const RealVector& etai  = m_collisions.etai(T, T, nd, p_x);
//        const VectorWrapper<double> x = asVector(p_x, ns);
//
//        double fac;
//
//        // Form the symmetric positive definite system matrix
//        int k = ns-nh;
//        for (int i = k; i < ns; ++i)
//            m_sys(i-k,i-k) = std::max(x(i) * x(i), 1.0E-99) / etai(i);
//
//        for (int i = k; i < ns; ++i) {
//            for (int j = i+1; j < ns; ++j) {
//                fac = std::max(x(i) * x(j), 1.0E-99) /
//                    (nDij(i,j) * (mi(i) + mi(j)));
//                m_sys(i-k,j-k) =  fac * (1.2 * Astar(i,j) - 2.0);
//                m_sys(i-k,i-k) += fac * (1.2 * mi(j) / mi(i) * Astar(i,j) + 2.0);
//                m_sys(j-k,j-k) += fac * (1.2 * mi(i) / mi(j) * Astar(i,j) + 2.0);
//            }
//        }
//
//        // Solve the linear system (this process is determined by the subclass)
//        static_cast<Implementation&>(*this).solve(x(k,ns));
//
//        // Finally compute the dot product of alpha and X
//        return dot(x(k,ns), m_alpha);
//    }
//
//protected:
//
//    RealSymMat m_sys;
//    RealVector m_alpha;
//
//}; // ViscosityCG
//
//
///**
// * Implements the conjugate gradient algorithm for solving the viscosity
// * transport system.
// */
//class ViscosityCG : public ViscositySPD<ViscosityCG>
//{
//    using ViscositySPD<ViscosityCG>::m_alpha;
//    using ViscositySPD<ViscosityCG>::m_sys;
//
//public:
//
//    ViscosityCG(CollisionDB& collisions)
//        : ViscositySPD<ViscosityCG>::ViscositySPD(collisions)
//    { }
//
//    template <typename E>
//    void solve(const VecExpr<double, E>& x)
//    {
//        DiagonalPreconditioner<double> M(m_sys);
//        cg(m_sys, m_alpha, x, M, 5);
//    }
//};
//
//// Register this algorithm with the other ViscosityAlgorithm types
//Config::ObjectProvider<ViscosityCG, ViscosityAlgorithm> viscosityCG("CG");
//
//
///**
// * Implements the LDL^T solution method for solving the viscosity transport
// * system.
// */
//class ViscosityLDLT : public ViscositySPD<ViscosityLDLT>
//{
//    using ViscositySPD<ViscosityLDLT>::m_alpha;
//    using ViscositySPD<ViscosityLDLT>::m_sys;
//
//public:
//
//    ViscosityLDLT(CollisionDB& collisions)
//        : ViscositySPD<ViscosityLDLT>::ViscositySPD(collisions)
//    { }
//
//    template <typename E>
//    void solve(const VecExpr<double, E>& x)
//    {
//        m_ldlt.setMatrix(m_sys);
//        m_ldlt.solve(m_alpha, x);
//    }
//
//private:
//
//    LDLT<double> m_ldlt;
//};
//
//// Register this algorithm with the other ViscosityAlgorithm types
//Config::ObjectProvider<ViscosityLDLT, ViscosityAlgorithm> viscosityLDLT("LDLT");


    } // namespace Transport
} // namespace Mutation

