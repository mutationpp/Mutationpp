/**
 * @file cg.h
 *
 * @brief Conjugate Gradient solver.
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

#ifndef NUMERICS_CG_H
#define NUMERICS_CG_H

#include "Preconditioner.h"

namespace Mutation {
    namespace Numerics {

/**
 * @ingroup itsol
 *
 * Solves the linear system \f$Ax = b\f$ using the Preconditioned Conjugate
 * Gradient method. \f$A\f$ is an NxN symmetric positive definite matrix and
 * \f$x\f$ and \f$b\f$ are Nx1 column vectors.
 */
template <typename T, typename E1, typename E2, typename Preconditioner>
void cg(const SymMatExpr<T, E1>& A, Vector<T>& x, const VecExpr<T, E2>& b, 
        const Preconditioner& M, const int max_iter = 100, 
        const double tolerance = 1.0e-10)
{
    T norm_b = b.norm2();
    T alpha, beta, rho0 = 0.0, rho1, resnorm;
    Vector<T> p(x.size()), q(x.size()), z(x.size());
    Vector<T> r = b - A * x;
    
    // Check for an already solved system
    if ((resnorm = r.norm2() / norm_b) <= tolerance) 
        return;
    
    // Perform CG iterations until convergence or reach max iterations
    for (int i = 1; i <= std::min(max_iter, static_cast<int>(A.rows())); ++i) {
        z = M.solve(r);
        rho1 = dot(r, z);
        
        if (i == 1) {
            p = z;
        } else {
            beta = rho1 / rho0;
            p = z + p * beta;
        }
        
        q = A * p;
        alpha = rho1 / dot(p, q);
        
        x = x + p * alpha;
        r = r - q * alpha;
        
        if ((resnorm = r.norm2() / norm_b) <= tolerance)     
            return;
       
        rho0 = rho1;
    } 
}

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_CG_H
