/**
 * @file ThermalConductivityChapmannEnskog.cpp
 *
 * @brief Provides Chapmann-Enskog thermal conductivity and thermal diffusion
 * ratios using LDLT and Conjugate-Gradient linear system solvers.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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
    ThermalConductivityChapmannEnskog(ThermalConductivityAlgorithm::ARGS args) :
        ThermalConductivityAlgorithm(args),
        m_sys(MatrixXd::Zero(args.nHeavy(), args.nHeavy())),
        m_x(args.nHeavy()),
        m_alpha(args.nHeavy())
    { }

    double thermalConductivity()
    {
        updateAlphas();
        return m_x.matrix().dot(m_alpha);
    }

    void thermalDiffusionRatios(double* const p_k)
    {
        // First solve the linear system for the alphas
        updateAlphas();

        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();
        const int k  = ns - nh;

        const ArrayXd& mi   = m_collisions.mass();
        const ArrayXd& nDij = m_collisions.nDij();
        const ArrayXd& Cst  = m_collisions.Cstij();

        // Compute heavy-particle thermal diffusion ratios
        double fac;
        m_sys.diagonal().setZero();
        for (int j = 0, si = 1; j < nh; ++j, ++si) {
            for (int i = j+1; i < nh; ++i, ++si) {
                fac = m_x(i)*m_x(j)/(mi(i+k)+mi(j+k))*(1.2*Cst(si)-1.)/nDij(si);
                m_sys(i,j) = fac*mi(i+k);
                m_sys(j,i) = fac*mi(j+k);
                m_sys(i,i) -= m_sys(j,i);
                m_sys(j,j) -= m_sys(i,j);
            }
        }

        p_k[0] = 0.0;
        Map<ArrayXd>(p_k+k,nh) = (m_sys*m_alpha)/KB;
    }

private:

    void updateAlphas()
    {
        const int ns = m_collisions.nSpecies();
        const int nh = m_collisions.nHeavy();

        // Collision data
        const ArrayXd& mi   = m_collisions.mass();
        const ArrayXd& Ast  = m_collisions.Astij();
        const ArrayXd& Bst  = m_collisions.Bstij();
        const ArrayXd& nDij = m_collisions.nDij();
        const ArrayXd& etai = m_collisions.etai();

        // Compute the system matrix
        double fac = 4.0 / (15.0 * KB);
        int k = ns - nh, ik, jk;
        m_x = Map<const ArrayXd>(m_collisions.thermo().X()+k,nh).max(1.0e-16);
        m_sys.diagonal().array() = fac*m_x*m_x*mi.tail(nh)/etai;

        double miij;
        double mjij;
        for (int j = 0, index = 1; j < nh; ++j, ++index) {
            jk = j+k;
            for (int i = j+1; i < nh; ++i, ++index) {
                ik = i+k;
                miij = mi(ik) / (mi(ik) + mi(jk));
                mjij = mi(jk) / (mi(ik) + mi(jk));
                fac  = m_x(i) * m_x(j) / (nDij(index) * 25.0 * KB);
                m_sys(i,j) = fac * miij * mjij *
                    (16.0 * Ast(index) + 12.0 * Bst(index) - 55.0);
                m_sys(i,i) += fac * (miij * (30.0 * miij + 16.0 * mjij *
                    Ast(index)) + mjij * mjij * (25.0 - 12.0 * Bst(index)));
                m_sys(j,j) += fac * (mjij * (30.0 * mjij + 16.0 * miij *
                    Ast(index)) + miij * miij * (25.0 - 12.0 * Bst(index)));
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

