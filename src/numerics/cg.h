#ifndef NUMERICS_CG_H
#define NUMERICS_CG_H

#include "Preconditioner.h"

namespace Numerics {

/**
 * Solves the linear system Ax = b using the Preconditioned Conjugate Gradient 
 * method.
 *
 * A is an NxN symmetric positive definite matrix and x and b are Nx1 column
 * vectors.
 */
template <typename T, typename E1, typename E2, typename Preconditioner>
void cg(const SymMatExpr<T, E1>& A, Vector<T>& x, const VecExpr<T, E2>& b, 
        const Preconditioner& M, const int max_iter = 100, 
        const double tolerance = 1.0e-6)
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

#endif // NUMERICS_CG_H
