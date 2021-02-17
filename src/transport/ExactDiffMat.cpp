/**
 * @file ExcactDiffMat.cpp
 *
 * @brief Implements ExcactDiffMat class.
 */

/*
 * Copyright 2016-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "AutoRegistration.h"
#include "CollisionDB.h"
#include "DiffusionMatrix.h"

#include <Eigen/Dense>

namespace Mutation {
    namespace Transport {

/**
 * Exact diffusion matrix implementation.
 */
class ExcactDiffMat : public DiffusionMatrix
{
public:

    ExcactDiffMat(DiffusionMatrix::ARGS collisions)
        : DiffusionMatrix(collisions)
    { }

    /**
     * Computes the exact multicomponent diffusion matrix.
     */
    const Eigen::MatrixXd& diffusionMatrix()
    {
        const int ns = m_collisions.nSpecies();
        const int k  = ns - m_collisions.nHeavy();
        const double nd = m_collisions.thermo().numberDensity();

        //Eigen::Map<const Eigen::ArrayXd> X = m_collisions.X();
        //Eigen::Map<const Eigen::ArrayXd> Y = m_collisions.Y();

        static Eigen::ArrayXd X; X = m_collisions.X()+1.0e-16; X /= X.sum();
        static Eigen::ArrayXd Y; Y.resize(ns);
        m_collisions.thermo().convert<Thermodynamics::X_TO_Y>(X.data(), Y.data());


        // Compute singular system matrix (lower triangular part)
        double fac;
        m_Dij.triangularView<Eigen::Lower>() = Eigen::MatrixXd::Zero(ns,ns);

        // Electron-Heavy subsystem
        if (k == 1) {
            const Eigen::ArrayXd& nDei = m_collisions.nDei();
            for (int i = 1; i < ns; ++i) {
                fac = X(0)*X(i)/nDei(i)*nd;
                m_Dij(0,0) += fac;
                m_Dij(i,i) += fac;
                m_Dij(i,0) = -fac;

            }
        }

        // Heavy particle subsystem
        const Eigen::MatrixXd& nDij = m_collisions.nDij();
        for (int j = k, si = 1; j < ns; ++j, ++si) {
            for (int i = j+1; i < ns; ++i, ++si) {
                fac = X(i)*X(j)/nDij(si)*nd;
                m_Dij(i,i) += fac;
                m_Dij(j,j) += fac;
                m_Dij(i,j) = -fac;
            }
        }

        // Add mass balance relation to make matrix nonsigular
        m_Dij.selfadjointView<Eigen::Lower>().rankUpdate(
            Y.matrix(), nd/nDij.diagonal().mean());

        static Eigen::LDLT<Eigen::MatrixXd, Eigen::Lower> ldlt;
        ldlt.compute(m_Dij.bottomRightCorner(ns-k,ns-k));

        static Eigen::VectorXd alpha; alpha.resize(ns-k);
        static Eigen::VectorXd b; b.resize(ns-k);
        b.array() = Y.tail(ns-k);

        for (int i = k; i < ns ; ++i ){
            b(i-k) += 1.0;
            alpha = ldlt.solve(b);
            b(i-k) -= 1.0;

            for (int j = i; j < ns; j++ ){
                m_Dij(i,j) = alpha(j-k);
                m_Dij(j,i) = alpha(j-k);
            }
        }

        return m_Dij;
    }

}; // ExcactDiffMat

// Register this algorithm
Mutation::Utilities::Config::ObjectProvider<
    ExcactDiffMat, DiffusionMatrix> exact_dm("Exact");

    } // namespace Transport
} // namespace Mutation

