/**
 * @file RamshawDiffMat.cpp
 *
 * @brief Implements RamshawDiffMat class.
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
 * Ramshaw diffusion matrix implementation.
 */
class RamshawDiffMat : public DiffusionMatrix
{
public:

    RamshawDiffMat(DiffusionMatrix::ARGS collisions)
        : DiffusionMatrix(collisions)
    { }

    /**
     * Computes the multicomponent diffusion coefficient matrix using the Fick
     * approximation with Ramshaw's correction,
     * \f[
     * D_{ij} = \frac{(\delta_{ij}-y_j)}{x_j}\frac{(1-y_j)}{(1-x_j)}D_{jm},
     * \f]
     * where \f$D_{jm}\f$ is average diffusion coefficient.
     * @see Mixture::averageDiffusionCoeffs().
     *
     * @todo This method currently breaks the symmetry property of the diffusion
     * matrix.  It should be reformulated to return a symmetric matrix which
     * would be faster to compute and use.
     */
    const Eigen::MatrixXd& diffusionMatrix()
    {
        const int ns = m_collisions.nSpecies();

        Eigen::Map<const Eigen::ArrayXd> X = m_collisions.X();
        Eigen::Map<const Eigen::ArrayXd> Y = m_collisions.Y();

        // Compute average diffusion coefficients without 1-X(j) term
        const Eigen::ArrayXd& Dim = m_collisions.Dim(false);

        // Form the matrix
        for (int j = 0; j < ns; ++j) {
            m_Dij.col(j).fill(-Y[j]/X[j]*(1.-Y[j])*Dim(j));
            m_Dij(j,j) -= m_Dij(j,j)/Y[j];
        }

        return m_Dij;
    }

}; // RamshawDiffMat

// Register this algorithm
Mutation::Utilities::Config::ObjectProvider<
    RamshawDiffMat, DiffusionMatrix> ramshaw_dm("Ramshaw");

    } // namespace Transport
} // namespace Mutation

