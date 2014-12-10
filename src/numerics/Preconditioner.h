/**
 * @file Preconditioner.h
 *
 * @brief Provides various preconditioners to be used in the linear solvers.
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

#ifndef NUMERICS_PRECONDITIONER_H
#define NUMERICS_PRECONDITIONER_H

#include "Vector.h"
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {


template <typename T>
class DiagonalPreconditioner
{
public:
    template <typename E>
    DiagonalPreconditioner(const MatExpr<T, E>& mat)
        : m_vec(mat.rows())
    {
        for (size_t i = 0; i < m_vec.size(); ++i)
            m_vec(i) = (mat(i,i) != 0.0 ? 1.0 / mat(i,i) : 1.0);
    }
    
    template <typename E>
    inline const VecMultVec<T, Vector<T>, E> solve(const VecExpr<T, E>& rhs) const {
        return m_vec * rhs;
    }
private:
    Vector<T> m_vec;   
};

template <typename T>
class IdentityPreconditioner
{
public:
    IdentityPreconditioner() { }
    
    template <typename E>
    inline const VecExpr<T,E>& solve(const VecExpr<T, E>& rhs) const {
        return rhs;
    }
};

    } // namespace Numerics
} // namespace Mutation


#endif // NUMERICS_PRECONDITIONER_H
