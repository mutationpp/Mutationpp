#ifndef NUMERICS_GMRES_H
#define NUMERICS_GMRES_H

#include <utility>
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {

/**
 * Solves the linear system MAx = b using the restarted, preconditioned
 * Generalized Minimum Residual (GMRES(m)) method.
 */
template <typename T, typename E1, typename E2, typename Preconditioner>
std::pair<int, double> gmres(
    const MatExpr<T, E1>& A, Vector<T>& x, const VecExpr<T, E2>& b, 
    const Preconditioner& M, const double tol = 1.0e-8, const int restart = 10, 
    const int max_iter = 100)
{
    unsigned int iter = 0;
    T bnrm2 = b.norm2();
    if (bnrm2 == static_cast<T>(0))
        bnrm2 = static_cast<T>(1);
        
    Vector<T> r = M.solve(b);
    T error = r.norm2()/bnrm2, temp;
    x = static_cast<T>(0);
    
    if (error < tol)
        return std::make_pair(0, error);
    
    const int n = A.rows();
    const int m = std::min(n, restart);
    
    std::vector< Vector<T> > V(m+1);
    
    Matrix<T> H(m+1,m);
    Vector<T> cs(m), sn(m), e1(m+1), s, w;
    e1(0) = 1.0;
    
    // Begin iteration
    while (iter < std::min(n, max_iter)) {
        r = M.solve(b - A*x);
        T rnorm = r.norm2();
        V[0] = r / rnorm;
        s = e1*rnorm;
        
        for (int i = 0; i < m; ++i, ++iter) {
            w = M.solve(A*V[i]);
            
            for (int k = 0; k <= i; ++k) {
                H(k,i) = dot(w, V[k]);
                w -= H(k,i)*V[k];
            }
            
            H(i+1,i) = w.norm2();
            V[i+1] = w / H(i+1,i);
            
            // Apply Givens rotations
            for (int k = 0; k < i; ++k) {
                temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
                H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
                H(k,i)   = temp;
            }
            
            // Compute the Given rotation for next column
            if (H(i+1,i) == static_cast<T>(0)) {
                cs(i) = static_cast<T>(1);
                sn(i) = static_cast<T>(0);
            } else if (std::abs(H(i+1,i)) > std::abs(H(i,i))) {
                temp = H(i,i) / H(i+1,i);
                sn(i) = 1.0/std::sqrt(1.0+temp*temp);
                cs(i) = temp*sn(i);
            } else {
                temp = H(i+1,i) / H(i,i);
                cs(i) = 1.0/std::sqrt(1.0+temp*temp);
                sn(i) = temp*cs(i);
            }
            
            s(i+1) = -sn(i)*s(i);
            s(i) = cs(i)*s(i);
            H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
            H(i+1,i) = static_cast<T>(0);
            
            error = std::abs(s(i+1)) / bnrm2;
            if (error <= tol) {
                update(i+1, H, s, V, x);
                return std::make_pair(iter, error);
            }
        }
        
        update(m, H, s, V, x);
        r = M.solve(b - A*x);
        s(m) = r.norm2();
        
        error = s(m) / bnrm2;
        if (error < tol)
            break;
    }
    
    return std::make_pair(iter, error);
}

/**
 * Helper function for the gmres function above.
 */
template <typename T>
void update(const int n, const Matrix<T>& H, const Vector<T>& s, 
    const std::vector< Vector<T> >& V, Vector<T>& x)
{
    // Back substitution to solve Hy = s
    Vector<T> y = s(0,n);
    y(n-1) /= H(n-1,n-1);
    for (int i = n-2; i >=0; --i) {
        for (int k = i+1; k < n; ++k)
            y(i) -= H(i,k)*y(k);
        y(i) /= H(i,i);
    }
    
    // Update x += Vy;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < V[0].size(); ++i)
            x(i) += V[j](i)*y(j);
}

    }
}

#endif // NUMERICS_GMRES_H

