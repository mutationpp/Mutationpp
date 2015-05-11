/**
 * @file LDLT.h
 *
 * @brief Provides LDLT class.
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

#ifndef NUMERICS_LDLT_H
#define NUMERICS_LDLT_H

#include "Vector.h"
#include "Matrix.h"

#include <iostream>
#include <cstdlib>

namespace Mutation {
    namespace Numerics {


/**
 * @ingroup dirsol
 * Performs the LDLT decomposition of an NxN symmetric positive definite
 * matrix A = LDL' and provides the solution algorithm for solving the SPD 
 * system Ax = b.
 */
template <typename T>
class LDLT
{
public:

    /**
     * Empty constructor, note that if you try and call solve() without setting
     * the matrix decomposition first using setMatrix(), the result will be an
     * error.
     */
    LDLT() {};

    /**
     * Constructs the LDLT decomposition A = LDL'.
     */
    template <typename E>
    LDLT(const SymMatExpr<T, E>& mat);
    
    /**
     * Constructs a new LDLT decomposition using the new matrix.  Returns false
     * if the matrix is numerically singular.
     */
    template <typename E>
    bool setMatrix(const SymMatExpr<T, E>& mat);
    
    /**
     * Solves the SPD system Ax = b via forward and backward substitutions with
     * the LDLT decomposition.
     */
    template <typename E>
    void solve(Vector<T>& x, const VecExpr<T, E>& b) const;
    
private:

    SymmetricMatrix<T> m_L;
    Vector<T>          m_work;
    
};


template <typename T>
template <typename E>
LDLT<T>::LDLT(const SymMatExpr<T, E>& mat)
{
    setMatrix(mat);
}

template <typename T>
template <typename E>
bool LDLT<T>::setMatrix(const SymMatExpr<T, E>& mat)
{
    m_L = mat;
    const size_t n = m_L.rows();
    m_work.resize(n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j)
            m_work(j) = m_L(j,j) * m_L(i,j);
        
        double& Lii = m_L(i,i);
        for (int j = 0; j < i; ++j)
            Lii -= m_work(j) * m_L(i,j);
        
        if (Lii == static_cast<T>(0)) {
            //std::cerr << "Calling LDLT decomposition with singular matrix!" << std::endl;
            //exit(1);
        	return false;
        }
        
        for (int k = i+1; k < n; ++k) {
            double& Lki = m_L(k,i);
            for (int j = 0; j < i; ++j)
                Lki -= m_L(k,j) * m_work(j);
            Lki /= Lii;
        }
    }

    return true;
}

template <typename T>
template <typename E>
void LDLT<T>::solve(Vector<T>& x, const VecExpr<T, E>& b) const
{
    const size_t n = m_L.rows();
    
    // Forward substitution
    x = b;
    for (int i = 0; i < n; ++i) {
        T& xi = x(i);
        for (int j = i+1; j < n; ++j)
            x(j) -= xi * m_L(i,j);
    }
    
    // Invert D
    for (int i = 0; i < n; ++i)
        x(i) /= m_L(i,i);
    
    // Backward substitution
    for (int i = n-2; i >= 0; --i) {
        T& xi = x(i);
        for (int j = i+1; j < n; ++j)
            xi -= x(j) * m_L(i,j);
    }
}

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_LDLT_H
