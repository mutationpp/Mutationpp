/**
 * @file DiffusionMatrix.h
 *
 * @brief Provides DiffusionMatrix base class.
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

#ifndef DIFFUSION_MATRIX_H
#define DIFFUSION_MATRIX_H

#include <Eigen/Dense>
#include "CollisionDB.h"

namespace Mutation {
    namespace Transport {

/**
 * Abstract base class for all diffusion matrix algorithms which allows for self
 * registration of concrete types.
 */
class DiffusionMatrix
{
public:

    // Required for self registering diffusion matrix algorithms
    typedef CollisionDB& ARGS;

    /// Returns name of this type.
    static std::string typeName() { return "DiffusionMatrix"; }

    /// Constructor taking a CollisionDB object.
    DiffusionMatrix(ARGS collisions)
        : m_collisions(collisions),
          m_Dij(collisions.nSpecies(), collisions.nSpecies())
    { }

    /// Destructor.
    virtual ~DiffusionMatrix() { }

    /// Returns the multicomponent diffusion matrix.
    virtual const Eigen::MatrixXd& diffusionMatrix() = 0;

protected:

    CollisionDB& m_collisions;
    Eigen::MatrixXd m_Dij;

}; // class DiffusionMatrix

    } // namespace Transport
} // namespace Mutation

#endif // DIFFUSION_MATRIX_H
