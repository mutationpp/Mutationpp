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
        CollisionDB& collisions)
        : m_thermo(thermo), m_collisions(collisions), 
          m_Di(thermo.nSpecies()), m_D(thermo.nSpecies(), thermo.nSpecies())
    { }
    
    /**
     * Computes the multicomponent diffusion coefficient matrix.
     *
     * @todo This method currently breaks the symmetry property of the diffusion
     * matrix.  It should be reformulated to return a symmetric matrix which
     * would be faster to compute and use.
     */
    const Eigen::MatrixXd& diffusionMatrix(
        const double T, const double Te, const double nd, const double *const p_x)
    {
        const int ns = m_thermo.nSpecies();
    
        // First step is to compute nDij and Y
        const Eigen::MatrixXd& nDij =
            m_collisions.nDij(T, Te, nd, p_x);
        double* Y = new double [ns];
        m_thermo.convert<Mutation::Thermodynamics::X_TO_Y>(p_x, Y);
        
        // Now we can compute the mixture averaged diffusion coefficient for the
        // ith species
        m_Di = Eigen::VectorXd::Zero(ns);
        for (int i = 0, index = 1; i < ns; ++i, ++index) {
            for (int j = i + 1; j < ns; ++j, ++index) {
                m_Di(i) += p_x[j] / nDij(index);
                m_Di(j) += p_x[i] / nDij(index);
            }
        }
        
        for (int i = 0; i < ns; ++i)
            m_Di(i) = (1.0 - Y[i]) / (m_Di(i) * p_x[i] * nd);
        
        // Form the diffusion matrix
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                m_D(i,j) = -Y[j] * m_Di(j);
            }
            m_D(i,i) += m_Di(i);
        }
        
        delete [] Y;
        return m_D;
    }
    
private:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    CollisionDB& m_collisions;
    Eigen::VectorXd m_Di;
    Eigen::MatrixXd m_D;
        
}; // class Ramshaw

    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_RANSHAW_H

