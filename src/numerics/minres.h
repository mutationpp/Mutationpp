#ifndef NUMERICS_MINRES_H
#define NUMERICS_MINRES_H

#include <utility>
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {

/**
 * Solves the symmetric indefinite (possibly singular) system Ax = b for x using
 * the MINRES algorithm.  This implementation is an adaptation of the MATLAB 
 * script provided in the following reference.
 *
 * "Polynomial based iteration methods for symmetric linear systems", Bernd 
 * Fischer, Wiley & Teubner, New York, Leipzig, 1996.
 *
 * @param A  the symmetric matrix
 * @param x  the solution vector (on return)
 * @param b  the right-hand-side vector
 * @param max_iter  maximum number of iterations to take
 * @param tol  the relative residual tolerance at which to stop iterating
 *
 * @return ret.first is the number of iterations used and ret.second is final
 * relative residual (may not be less than tol if MINRES fails to converge)
 */
template <typename T, typename E1, typename E2>
std::pair<int, double> minres(
    const SymMatExpr<T, E1>& A, Vector<T>& x, const VecExpr<T, E2>& b,
    const int max_iter = 100, const double tol = 1.0e-6)
{
    assert(A.rows() == A.cols());
    assert(A.cols() == x.size());
    assert(A.rows() == b.size());

    int iter = 0;
    int n = x.size();
    
    Vector<T> v(n); 
    Vector<T> vhat = b;
    Vector<T> vold = v;
    T beta = vhat.norm2(), betaold;
    T c    = static_cast<T>(1);
    T cold = static_cast<T>(1), coold;
    T s    = static_cast<T>(0);
    T sold = static_cast<T>(0), soold;
    T rhat, r1, r2, r3, norm = beta, norm0 = beta, alpha;
    
    Vector<T> w(n);
    Vector<T> wold = w, woold(n);
    T eta = beta;
    x = static_cast<T>(0);
    
    while (iter < std::min(n, max_iter) && (norm / norm0 > tol)) {
        iter++;
        
        // Lanczos
        vold = v;  v = vhat / beta;
        vhat = A*v;  alpha = dot(v, vhat);  vhat -= alpha*v + beta*vold;
        betaold = beta;  beta = vhat.norm2();
        
        // QR factorization
        coold = cold; cold = c; soold = sold; sold = s;
        rhat = cold*alpha - coold*sold*betaold; 
        r1   = std::sqrt(rhat*rhat + beta*beta);
        r2   = sold*alpha + coold*cold*betaold;
        r3   = soold*betaold;
        
        // Givens rotation
        c = rhat / r1;
        s = beta / r1;
        
        // Update
        woold = wold; wold = w;
        w = (v - r3*woold - r2*wold) / r1;
        x += c*eta*w; norm *= std::abs(s);
        eta *= -s;
    }
    
    return std::make_pair(iter, norm / norm0);    
}

    } // Numerics
} // Mutation

#endif // NUMERICS_MINRES_H


