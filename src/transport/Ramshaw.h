/**
 * @file Ramshaw.h
 *
 * @brief Provides Ramshaw class.
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

#ifndef TRANSPORT_RANSHAW_H
#define TRANSPORT_RANSHAW_H

#include "CollisionDB.h"
#include "Thermodynamics.h"

#include "Eigen/Dense"

namespace Mutation {
    namespace Transport {

class Ramshaw
{
public:
    
    Ramshaw(
        const Mutation::Thermodynamics::Thermodynamics& thermo, 
        CollisionDB& collisions) :
        m_thermo(thermo), m_collisions(collisions),
        m_D(thermo.nSpecies(), thermo.nSpecies())
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
    const Eigen::MatrixXd& diffusionMatrix(
        const double T, const double Te, const double nd, const double *const p_x)
    {
        const int ns = m_thermo.nSpecies();
    
        // First step is to compute X and Y with tolerance on X
        static Eigen::ArrayXd X(ns);
        static Eigen::ArrayXd Y(ns);
        //X = Eigen::Map<const Eigen::ArrayXd>(p_x, ns)+1.0e-16;
        X = Eigen::Map<const Eigen::ArrayXd>(p_x, ns);//.max(1.0e-30);
        X /= X.sum();
        m_thermo.convert<Mutation::Thermodynamics::X_TO_Y>(X.data(), Y.data());
        
        // Compute average diffusion coefficients
        const Eigen::ArrayXd& Dim = m_collisions.Dim(T, Te, nd, p_x);

        // Form the matrix
        for (int j = 0; j < ns; ++j) {
            m_D.col(j).fill(-Y[j]/X[j]*(1.-Y[j])/(1.-X[j])*Dim(j));
            m_D(j,j) -= m_D(j,j)/Y[j];
        }

        return m_D;
    }
    
private:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    CollisionDB& m_collisions;
    Eigen::MatrixXd m_D;
        
}; // class Ramshaw

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_RANSHAW_H

