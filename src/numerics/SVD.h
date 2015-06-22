/**
 * @file SVD.h
 *
 * @brief Implements the SVD class which computes the singular value
 * decomposition of a matrix.
 *
 * @see class SVD
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

#ifndef NUMERICS_SVD_H
#define NUMERICS_SVD_H

#include <cassert>
#include <cmath>

#include "NumConst.h"
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {

/**
 * @ingroup dirsol
 * Represents the singular value decomposition of a real MxN matrix A = U*S*V',
 * where U is a real MxN orthogonal matrix and V is a real NxN orthogonal matrix
 * and S is the diagonal matrix of singular values listed in decreasing order.
 *
 * <p>
 * (This class is an adaptation of the JAMA/C++ SVD class developed jointly by 
 * Mathworks and NIST.  See http://math.nist.gov/javanumerics/jama.)
 *
 * @todo The actual SVD algorithm was predominantly taken from the JAMA/C++ code
 * and is therefore not necessarily completely portable in terms of efficiency 
 * for the Matrix and Vector clases used here.  Therefore, care should be taken
 * to inspect the SVD constructor to make sure that matrix and vector 
 * multiplies, etc., are handled in the appropriate order (ie: row major 
 * ordering).  From an initial glance it looks like they are not.  In addition,
 * some optimizations are made based on column ordering (ie: copying the row of
 * A into e during bidiagonalization) which should be either undone or reworked
 * for the appropriate cases.
 *
 * @todo The solve(), U(), and V() commands all currently assume that U and V
 * are already computed.  The real value of putting the SVD into a class is that
 * it could delay the computation of U and V until the user requests them to be
 * computed.  In solve, instead of using U directly, U'b could be computed using
 * the factored Householder vectors if U is not computed already.  These 
 * changes could drastically improve performance for some problems but are not
 * necessary unless the SVD becomes a bottleneck during performance testing.
 */
template <typename Real>
class SVD
{
public:
    
    /**
     * Computes the singular value decomposition of A = USV' in factored form.
     */
    template <typename E>
    SVD(const MatExpr<Real, E>& A,
        const bool wantu = true, 
        const bool wantv = true);
    
    /**
     * Returns the singular values of A.
     */
    const Vector<Real> &singularValues() const {
        return m_s;
    }
    
    /**
     * Returns the left orthogonal matrix U.
     */
    const Matrix<Real> &U() const {
        return m_U;
    }
    
    /**
     * Returns the right orthogonal matrix V.
     */
    const Matrix<Real> &V() const {
        return m_V;
    }
    
    /**
     * Returns the rank of A.
     */
    int rank() const {
        return m_rank;
    }
    
    /**
     * Returns the two-norm of A.
     */
    Real norm() const {
        return m_s(0);
    }
    
    /**
     * Returns the two-norm condition number of A.
     */
    Real cond() const {
        return m_s(0) / m_s(min(m_A.rows(),m_A.cols()));
    }
    
    /**
     * Determines the minimum two norm solution of the least squares problem
     * A*x = b.
     */
    void solve(Vector<Real> &x, Vector<Real> &b) const;
    
    /**
     * Determines the minimum two norm solution of the least squares problem
     * A*x = b.
     */
    void solve(Real* const x, Real* const b) const;
    
    /**
     * Solves the system A'Ax = b via VS^2V'x = b
     */
    void solveATA(Vector<Real> &x, Vector<Real> &b);

private:

    int m_rank;

    Matrix<Real> m_A;
    Vector<Real> m_s;
    Matrix<Real> m_U;
    Matrix<Real> m_V;
};

template <typename Real>
Real hypot(const Real &a, const Real &b)
{
    if (a == NumConst<Real>::zero)
        return std::abs(b);
    
    if (std::abs(a) > std::abs(b)) {
        Real c = b / a;
        return std::abs(a) * std::sqrt(NumConst<Real>::one + c * c);
    } else {
        Real c = a / b;
        return std::abs(b) * std::sqrt(NumConst<Real>::one + c * c);
    }
}

template <typename Real>
template <typename E>
SVD<Real>::SVD(const MatExpr<Real, E> &A, const bool wantu, const bool wantv)
    : m_A(A)
{
    // Define size constants
    const int m   = A.rows();           // number of rows in A
    const int n   = A.cols();           // number of columns in A
    //const int nu  = min(m,n);           // number of columns in U
    const int nct = std::min(m-1,n);         // number of column x-forms for bidiag
    const int nrt = std::max(0, std::min(n-2,m)); // number of row x-forms for bidiag
        
    m_s = Vector<Real>(std::min(m+1,n));
    m_U = Matrix<Real>(m, m);  //**** 65
    m_V = Matrix<Real>(n, n);

    Vector<Real> e(n);
    Vector<Real> work(m);
	
	int i, j, k;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.
    for (k = 0; k < std::max(nct,nrt); k++) {
        if (k < nct) {
            // Compute the transformation for the k-th column and
            // place the k-th diagonal in m_s(k).
            // Compute 2-norm of k-th column without under/overflow.
            m_s(k) = NumConst<Real>::zero;
            for (i = k; i < m; i++) {
                m_s(k) = hypot(m_s(k), m_A(i,k));
            }
            if (m_s(k) != NumConst<Real>::zero) {
                if (m_A(k,k) < NumConst<Real>::zero) {
                    m_s(k) = -m_s(k);
                }
                for (i = k; i < m; i++) {
                    m_A(i,k) /= m_s(k);
                }
                m_A(k,k) += NumConst<Real>::one;
            }
            m_s(k) = -m_s(k);
        }
        for (j = k+1; j < n; j++) {
            if ((k < nct) && (m_s(k) != NumConst<Real>::zero))  {
                // Apply the transformation.
                double t = 0;
                for (i = k; i < m; i++) {
                    t += m_A(i,k) * m_A(i,j);
                }
                t = -t/m_A(k,k);
                for (i = k; i < m; i++) {
                    m_A(i,j) += t * m_A(i,k);
                }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.
            e(j) = m_A(k,j);
        }
        if (wantu & (k < nct)) {
            // Place the transformation in U for subsequent back
            // multiplication.
            for (i = k; i < m; i++) {
                m_U(i,k) = m_A(i,k);
            }
        }
        if (k < nrt) {
            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e(k).
            // Compute 2-norm without under/overflow.
            e(k) = NumConst<Real>::zero;
            for (i = k+1; i < n; i++) {
                e(k) = hypot(e(k),e(i));
            }
            if (e(k) != NumConst<Real>::zero) {
                if (e(k+1) < NumConst<Real>::zero) {
                    e(k) = -e(k);
                }
                for (i = k+1; i < n; i++) {
                    e(i) /= e(k);
                }
                e(k+1) += NumConst<Real>::one;
            }
            e(k) = -e(k);
            
            if ((k+1 < m) & (e(k) != NumConst<Real>::zero)) {
                // Apply the transformation.
                for (i = k+1; i < m; i++) {
                    work(i) = NumConst<Real>::zero;
                }
                for (j = k+1; j < n; j++) {
                    for (i = k+1; i < m; i++) {
                        work(i) += e(j) * m_A(i,j);
                    }
                }
                for (j = k+1; j < n; j++) {
                    double t = -e(j) / e(k+1);
                    for (i = k+1; i < m; i++) {
                        m_A(i,j) += t * work(i);
                    }
                }
            }
            if (wantv) {
                // Place the transformation in V for subsequent
                // back multiplication.
                for (i = k+1; i < n; i++) {
                    m_V(i,k) = e(i);
                }
            }
        }
    }

    // Set up the final bidiagonal matrix of order p.
    int p = std::min(n,m+1);
    if (nct < n)   m_s(nct) = m_A(nct,nct);        
    if (m < p)     m_s(p-1) = NumConst<Real>::zero;
    if (nrt+1 < p) e(nrt)   = m_A(nrt,p-1);
    e(p-1) = NumConst<Real>::zero;

    // If required, generate U.
    if (wantu) {
        for (j = nct; j < m; j++) {
            for (i = 0; i < m; i++) {
                m_U(i,j) = NumConst<Real>::zero;
            }
            m_U(j,j) = NumConst<Real>::one;
        }
        for (k = nct-1; k >= 0; k--) {
            if (m_s(k) != NumConst<Real>::zero) {
                for (j = k+1; j < m; j++) {
                    double t = NumConst<Real>::zero;
                    for (i = k; i < m; i++) {
                        t += m_U(i,k) * m_U(i,j);
                    }
                    t = -t/m_U(k,k);
                    for (i = k; i < m; i++) {
                        m_U(i,j) += t * m_U(i,k);
                    }
                }
                for (i = k; i < m; i++ ) {
                    m_U(i,k) = -m_U(i,k);
                }
                m_U(k,k) = NumConst<Real>::one + m_U(k,k);
                for (i = 0; i < k-1; i++) {
                    m_U(i,k) = NumConst<Real>::zero;
                }
            } else {
                for (i = 0; i < m; i++) {
                    m_U(i,k) = NumConst<Real>::zero;
                }
                m_U(k,k) = NumConst<Real>::one;
            }
        }
    }

    // If required, generate V.
    if (wantv) {
        for (k = n-1; k >= 0; k--) {
            if ((k < nrt) & (e(k) != NumConst<Real>::zero)) {
                for (j = k+1; j < n; j++) {
                    double t = NumConst<Real>::zero;
                    for (i = k+1; i < n; i++) {
                        t += m_V(i,k) * m_V(i,j);
                    }
                    t = -t / m_V(k+1,k);
                    for (i = k+1; i < n; i++) {
                        m_V(i,j) += t * m_V(i,k);
                    }
                }
            }
            for (i = 0; i < n; i++) {
                m_V(i,k) = NumConst<Real>::zero;
            }
            m_V(k,k) = NumConst<Real>::one;
        }
    }

    // Main iteration loop for the singular values.
    int pp = p-1;
    int iter = 0;
    
    while (p > 0) {
		int kase = 0;
		Real eps = NumConst<Real>::eps;
		Real tiny = eps*eps;

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e(k-1) are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e(k-1) is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p-2; k >= -1; k--) {
            if (k == -1) 
                break;
            
            if (std::abs(e(k)) <= tiny + eps*(std::abs(m_s(k)) + std::abs(m_s(k+1)))) {
               e(k) = NumConst<Real>::zero;
               break;
            }
        }
        if (k == p-2) {
            kase = 4;
        } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
                if (ks == k)
                    break;
                
                double t = (ks != p ? std::abs(e(ks)) : NumConst<Real>::zero) + 
                           (ks != k+1 ? std::abs(e(ks-1)) : NumConst<Real>::zero);
               
                if (std::abs(m_s(ks)) <= tiny + eps * t)  {
                    m_s(ks) = NumConst<Real>::zero;
                    break;
                }
            }
            
            if (ks == k) {
                kase = 3;
            } else if (ks == p-1) {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }
        k++;

        // Perform the task indicated by kase.
        switch (kase) {
            // Deflate negligible s(p).
            case 1: {
                double f = e(p-2);
                e(p-2) = NumConst<Real>::zero;
                for (j = p-2; j >= k; j--) {
                    double t = hypot(m_s(j),f);
                    double cs = m_s(j) / t;
                    double sn = f / t;
                    m_s(j) = t;
                    
                    if (j != k) {
                        f = -sn * e(j-1);
                        e(j-1) = cs * e(j-1);
                    }
                    
                    if (wantv) {
                        for (i = 0; i < n; i++) {
                            t = cs * m_V(i,j) + sn * m_V(i,p-1);
                            m_V(i,p-1) = -sn * m_V(i,j) + cs * m_V(i,p-1);
                            m_V(i,j) = t;
                        }
                    }
                }
            }
            break;

            // Split at negligible s(k).
            case 2: {
                double f = e(k-1);
                e(k-1) = NumConst<Real>::zero;
                
                for (j = k; j < p; j++) {
                    double t = hypot(m_s(j), f);
                    double cs = m_s(j) / t;
                    double sn = f / t;
                    m_s(j) = t;
                    f = -sn * e(j);
                    e(j) = cs * e(j);
                  
                    if (wantu) {
                        for (i = 0; i < m; i++) {
                            t = cs * m_U(i,j) + sn * m_U(i,k-1);
                            m_U(i,k-1) = -sn * m_U(i,j) + cs * m_U(i,k-1);
                            m_U(i,j) = t;
                        }
                    }
                }
            }
            break;

            // Perform one qr step.
            case 3: {
                // Calculate the shift.   
                Real scale = std::max(std::max(std::max(std::max(
                               std::abs(m_s(p-1)),std::abs(m_s(p-2))),std::abs(e(p-2))), 
                               std::abs(m_s(k))),std::abs(e(k)));
                Real sp = m_s(p-1) / scale;
                Real spm1 = m_s(p-2) / scale;
                Real epm1 = e(p-2) / scale;
                Real sk = m_s(k) / scale;
                Real ek = e(k) / scale;
                Real b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) /
                           NumConst<Real>::two;
                Real c = (sp * epm1) * (sp * epm1);
                Real shift = NumConst<Real>::zero;
                
                if ((b != NumConst<Real>::zero) || 
                    (c != NumConst<Real>::zero)) {
                    shift = std::sqrt(b*b + c);
                    if (b < NumConst<Real>::zero) {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
               
                Real f = (sk + sp)*(sk - sp) + shift;
                Real g = sk * ek;
   
                // Chase zeros.   
                for (j = k; j < p-1; j++) {
                    Real t = hypot(f,g);
                    Real cs = f / t;
                    Real sn = g / t;
                  
                    if (j != k)
                        e(j-1) = t;

                    f = cs * m_s(j) + sn * e(j);
                    e(j) = cs * e(j) - sn * m_s(j);
                    g = sn * m_s(j+1);
                    m_s(j+1) = cs * m_s(j+1);
                    
                    if (wantv) {
                        for (i = 0; i < n; i++) {
                            t = cs * m_V(i,j) + sn * m_V(i,j+1);
                            m_V(i,j+1) = -sn * m_V(i,j) + cs * m_V(i,j+1);
                            m_V(i,j) = t;
                        }
                    }
                    
                    t = hypot(f,g);
                    cs = f / t;
                    sn = g / t;
                    m_s(j) = t;
                    f = cs * e(j) + sn * m_s(j+1);
                    m_s(j+1) = -sn * e(j) + cs * m_s(j+1);
                    g = sn * e(j+1);
                    e(j+1) = cs * e(j+1);
                    
                    if (wantu && (j < m-1)) {
                        for (i = 0; i < m; i++) {
                            t = cs * m_U(i,j) + sn * m_U(i,j+1);
                            m_U(i,j+1) = -sn * m_U(i,j) + cs * m_U(i,j+1);
                            m_U(i,j) = t;
                        }
                    }
                }
                e(p-2) = f;
                iter = iter + 1;
            }
            break;

            // Convergence.
            case 4: {
                // Make the singular values positive.   
                if (m_s(k) <= NumConst<Real>::zero) {
                    m_s(k) = (m_s(k) < NumConst<Real>::zero ? 
                                -m_s(k) : NumConst<Real>::zero);
                  
                    if (wantv) {
                        for (i = 0; i <= pp; i++) {
                            m_V(i,k) = -m_V(i,k);
                        }
                    }
                }
   
                // Order the singular values.   
                while (k < pp) {
                    if (m_s(k) >= m_s(k+1))
                        break;

                    double t = m_s(k);
                    m_s(k) = m_s(k+1);
                    m_s(k+1) = t;
                  
                    if (wantv && (k < n-1)) {
                        for (i = 0; i < n; i++) {
                            t = m_V(i,k+1); 
                            m_V(i,k+1) = m_V(i,k); 
                            m_V(i,k) = t;
                        }
                    }
                    
                    if (wantu && (k < m-1)) {
                        for (i = 0; i < m; i++) {
                            t = m_U(i,k+1); 
                            m_U(i,k+1) = m_U(i,k); 
                            m_U(i,k) = t;
                        }
                    }
                    
                    k++;
                }
                iter = 0;
                p--;
            }
            break;
        }
    }
    
    // Now determine and store the rank
    Real tol = m_s.size() * m_s(0) * NumConst<Real>::eps;
    m_rank = std::min(m,n);
        
    while (m_s(m_rank-1) < tol)
        m_rank--;
} // SVD()

template <typename Real>
void SVD<Real>::solve(Vector<Real> &x, Vector<Real> &b) const
{
    const int m = m_A.rows();
    const int n = m_A.cols();
    
    x(0,m_rank) = b(0,m) * m_U(0,m,0,m_rank);
    
    for (int k = 0; k < m_rank; ++k) {
        b(k) = x(k) / m_s(k);
        x(k) = NumConst<Real>::zero;
    }
    
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < m_rank; ++k)
            x(i) += m_V(i,k) * b(k);
}

template <typename Real>
void SVD<Real>::solve(Real* const x, Real* const b) const
{
    const int m = m_A.rows();
    const int n = m_A.cols();
    
    for (int i = 0; i < m_rank; ++i) {
        x[i] = 0.0;
        for (int j = 0; j < m; ++j)
            x[i] += b[j]*m_U(j,i);
    }
    
    for (int k = 0; k < m_rank; ++k) {
        b[k] = x[k] / m_s(k);
        x[k] = Real(0.0);
    }
    
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < m_rank; ++k)
            x[i] += m_V(i,k) * b[k];
}

template <typename Real>
void SVD<Real>::solveATA(Vector<Real> &x, Vector<Real> &b)
{
    x = b * m_V;
    b = x / square(m_s);
    x = m_V * b;
}


    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_SVD_H
