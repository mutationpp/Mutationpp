/**
 * @file ViscosityChapmannEnskog.cpp
 *
 * @brief Provides Chapmann-Enskog viscosity algorithms using LDLT and
 * Conjugate-Gradient linear system solvers.
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

#include "ViscosityAlgorithm.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

/**
 * Computes the viscosity using the full Chapmann-Enskog formulation.
 *
 * @param Solver an Eigen linear system solver
 */
template <template <typename, int> class Solver>
class ViscosityChapmannEnskog : public ViscosityAlgorithm
{
public:
	ViscosityChapmannEnskog(CollisionDB& collisions) :
		ViscosityAlgorithm(collisions),
		m_sys(MatrixXd::Zero(collisions.nHeavy(), collisions.nHeavy())),
		m_x(collisions.nHeavy()),
		m_alpha(collisions.nHeavy())
	{ }

	double viscosity(
	    double Th, double Te, double nd, const double *const p_x)
	{
		const int ns = m_collisions.nSpecies();
		const int nh = m_collisions.nHeavy();

		const ArrayXd&  mi    = m_collisions.mass();
		const ArrayXXd& Astar = m_collisions.Astar(Th, Te, nd, p_x);
		const ArrayXXd& nDij  = m_collisions.nDij(Th, Te, nd, p_x);
		const ArrayXd&  etai  = m_collisions.etai(Th, Te, nd, p_x);

		// Form the symmetric positive definite system matrix
		int k = ns-nh, ik, jk;
		double fac;

		m_x = Map<const ArrayXd>(p_x+k,nh).max(1.0e-30);
        m_sys.diagonal().array() = m_x * m_x / etai.tail(nh);

		for (int j = k; j < ns; ++j) {
			jk = j-k;
			for (int i = j+1; i < ns; ++i) {
				ik = i-k;
				fac = m_x(ik)*m_x(jk) / (nDij(i,j) * (mi(i) + mi(j)));
				m_sys(ik,jk) =  fac * (1.2 * Astar(i,j) - 2.0);
				m_sys(ik,ik) += fac * (1.2 * mi(j) / mi(i) * Astar(i,j) + 2.0);
				m_sys(jk,jk) += fac * (1.2 * mi(i) / mi(j) * Astar(i,j) + 2.0);
			}
		}

		// Solve the linear system
		m_solver.compute(m_sys);
		m_alpha = m_solver.solve(m_x.matrix());

		// Finally compute the dot product of alpha and X
		return m_x.matrix().dot(m_alpha);
	}

private:

	MatrixXd m_sys;
	ArrayXd m_x;
	VectorXd m_alpha;
	Solver<MatrixXd, Lower> m_solver;
};

// Register the Chapmann-Enskog solution using the LDLT decomposition
Config::ObjectProvider<
	ViscosityChapmannEnskog<LDLT>, ViscosityAlgorithm>
	visc_CE_LDLT("Chapmann-Enskog_LDLT");

// Register the Chapmann-Enskog solution using the Conjugate-Gradient algorithm.
// Note the use of a proxy class CG to reduce the template arguments to 2
template <typename MatrixType, int UpLo>
class CG : public ConjugateGradient<MatrixType, UpLo> { };
Config::ObjectProvider<
	ViscosityChapmannEnskog<CG>, ViscosityAlgorithm>
	visc_CE_CG("Chapmann-Enskog_CG");


    } // namespace Transport
} // namespace Mutation

