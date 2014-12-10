/**
 * @file gmres.h
 *
 * @brief Generalized Minimum Residual (GMRES) solver.
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

#ifndef NUMERICS_GMRES_H
#define NUMERICS_GMRES_H

#include <utility>
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {

/**
 * @ingroup itsol
 * Solves the linear system inv(M)Ax = inv(b) using the preconditioned
 * Generalized Minimum Residual (GMRES) method.
 */
template <typename T, typename E1, typename E2, typename P>
std::pair<int, double> gmres(
    const MatExpr<T, E1>& A, Vector<T>& x, const VecExpr<T, E2>& b,
    const P& M, const double tol = 1.0e-8, const int max_iter = -1)
{
    // Figure out sizing information
    const int n = A.rows();
    const int m = (max_iter < 0 ? n : max_iter);

    std::vector< Vector<T> > V(m+1);
    Matrix<T> H(m+1,m, 0.0);
    Vector<T> r(m+1, 0.0);
    Vector<T> c(m, 0.0);
    Vector<T> s(m, 0.0);

    // Initial guess is zero
    x = static_cast<T>(0.0);
    V[0] = M.solve(b);
    double norm0 = V[0].norm2();

    // First term of the Hessenberg system
    if (norm0 > static_cast<T>(0))
        V[0] /= norm0;
    else
        return std::make_pair(0, 0.0);

    // Iterative loop
    double norm = norm0, t;
    int i = 0, i1;
    r(0) = norm;

    while (i < m && norm > norm0 * tol) {
        i1 = i+1;
        V[i1] = M.solve(A*V[i]);

        // Modified Graham-Schmidt and Arnoldi step
        for (int j = 0; j < i1; ++j) {
            H(j,i) = dot(V[j], V[i1]);
            V[i1] -= H(j,i)*V[j];
        }

        H(i1,i) = V[i1].norm2();
        if (H(i1,i) != static_cast<T>(0))
            V[i1] /= H(i1,i);

        // Factorization of H
        for (int j = 0; j < i; ++j) {
            t = H(j,i);
            H(j,i)   =  c(j)*t + s(j)*H(j+1,i);
            H(j+1,i) = -s(j)*t + c(j)*H(j+1,i);
        }

        t = std::sqrt(H(i,i)*H(i,i)+H(i1,i)*H(i1,i));
        if (t == static_cast<T>(0.0))
            t = std::numeric_limits<T>::epsilon();

        // Next plane rotation
        c(i)  =  H(i,i)/t;
        s(i)  =  H(i1,i)/t;
        r(i1) = -s(i)*r(i);
        r(i)  =  c(i)*r(i);

        // Residual norm
        H(i,i) = c(i)*H(i,i) + s(i)*H(i1,i);
        H(i1,i) = static_cast<T>(0.0);
        norm = std::abs(r(i1));
        i = i1;
    }

    // Solve for upper triangular system
    for (int k = i-1; k >= 0; --k) {
        t = static_cast<T>(0);
        for (int j = k+1; j < i; ++j)
            t += H(k,j)*r(j);
        r(k) = (r(k) - t)/H(k,k);
    }

    // Form linear combination of V(:,i)'s to get solution
    for (int j = 0; j < i; ++j)
        x += r(j)*V[j];
    
    return std::make_pair(i,norm/norm0);
}

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_GMRES_H

