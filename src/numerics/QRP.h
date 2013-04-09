/**
 * @file QRP.h
 *
 * Implements the QRP class which performs the QR decomposition with column 
 * pivoting.
 *
 * @see class QR
 * @see class QRP
 *
 * @author J.B. Scoggins (jbscoggi@gmail.com)
 * @date   February 2, 2012
 */
 
#ifndef NUMERICS_QRP_H
#define NUMERICS_QRP_H

#include "QR.h"

namespace Mutation {
    namespace Numerics {

/**
 * Computes the QR decomposition with column pivoting of an MxN matrix A where
 * M >= N using Householder reflections.  The decomposition of AP = QR is such
 * that Q is MxM orthogonal, R is MxN upper triangular, and P is an NxN
 * permutation matrix.
 */
template <typename Real>
class QRP : public QR<Real>
{
public:

    /**
     * Constructs the QR decomposition with column pivoting of A.
     */
    QRP(const Matrix<Real> &A);
    
    /**
     * Returns the estimated rank of matrix A as deduced from the QRP algorithm.
     */
    const int &rank() {
        return m_rank;
    }
    
    /**
     * Applies the permutation matrix P to a Vector x, ie: x = Px.  If transpose
     * is true, then this method applies the transpose of the permutation matrix
     * x = P'x.
     */
    template <typename T>
    void pivot(Vector<T> &x, const bool transpose = false) const;
    
    /**
     * Returns a Vector<int> whos i'th entry is the column in AP which was
     * the i'th column in A.
     */
    Vector<int> columnOrder() const;
    
    /**
     * Solves the linear system Ax = b via RP'x = Q'b.
     */
    void solve(Vector<Real> &x, Vector<Real> &b) const;
    
    /**
     * Solves the system A'Ax = b via PR'RP'x = b.
     */
    void solveATA(Vector<Real> &x, Vector<Real>&b) const;
        
    // Templated base class requires following directives to find protected
    // members
    using QR<Real>::m_A;
    using QR<Real>::m_betas;

private:
    
    std::vector<int> m_pivots;
    int m_rank;
};

//==============================================================================

template <typename Real>
QRP<Real>::QRP(const Matrix<Real> &A)
    : QR<Real>(A, true), m_rank(0)
{
    const size_t m = m_A.rows();
    const size_t n = m_A.cols();
    
    assert( m >= n );
    
    size_t i, j;
    
    Vector<Real> v(n);
    Vector<Real> u(m);
    
    // Store the column norms in vector c
    Vector<Real> c(n);
    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            c(j) += m_A(i,j) * m_A(i,j);
    
    // Find column with largest norm
    IVP<Real> max = c.max();
    size_t k = max.index();
    Real tau = max.value();
    
    while (tau > static_cast<Real>(0)) {    
        m_pivots.push_back(k);
        
        // Perform swap if necessary
        if (m_rank != k) {            
            // Swap columns of A
            m_A.swapCols(m_rank, k);            
            // Swap elements of column norm vector
            c.swap(m_rank, k);
        }
        
        // Compute Householder vector and multiplier beta
        u(m_rank,m) = m_A(m_rank,m,m_rank,m_rank+1);
        house(u(m_rank,m), m_betas(m_rank));
        
        // Rank-1 update of A(r:m,r:n) = (I - beta*v*v')A(r:m,r:n)
        v(m_rank,n) = m_betas(m_rank) * (u(m_rank,m) * m_A(m_rank,m,m_rank,n));
        m_A(m_rank,m_rank) -= v(m_rank);
        
        if (m_rank + 1 < n)
            m_A(m_rank,m,m_rank+1,n) -= rank1(u(m_rank,m), v(m_rank+1,n));
        
        // Store Householder vector in zeroed out portion of A
        if (m_rank + 1 < m)
            m_A(m_rank+1,m,m_rank,m_rank+1) = u(m_rank+1,m);
        
        // Downdate column norm
        for (i = m_rank+1; i < n; ++i)
            c(i) -= m_A(m_rank,i) * m_A(m_rank,i);
        
        // Find the next pivot column
        if (m_rank < n-1) {
            max = c(m_rank+1,n).max();
            k   = max.index() + m_rank + 1;
            tau = max.value();
        } else {
            tau = static_cast<Real>(0);
        }
        
        // Increment rank    
        m_rank++;
    }
}

//==============================================================================

template <typename Real>
template <typename T>
void QRP<Real>::pivot(Vector<T> &x, const bool transpose) const
{
    const int npiv = m_pivots.size();
    
    if (transpose) {
        for(size_t i = 0; i < npiv; ++i)
            x.swap(i, m_pivots[i]);
    } else {
        for(int i = npiv-1; i >= 0; --i)
            x.swap(i, m_pivots[i]);
    }
}

//==============================================================================

template <typename Real>
Vector<int> QRP<Real>::columnOrder() const
{
    const int nc = m_A.cols();
    
    Vector<int> order(nc);
    for (int i = 0; i < nc; ++i)
        order(i) = i;
    
    pivot(order, true);    
    return order;
}

//==============================================================================

template <typename Real>
void QRP<Real>::solve(Vector<Real> &x, Vector<Real> &b) const
{
    // Update Q'b
    qTransposeB(b);
    
    // Back-substitution
    for (int i = m_A.cols()-1; i >= 0; --i) {
        for (int j = i+1; j < m_A.cols(); ++j)
            b(i) -= m_A(i,j) * x(j);
        x(i) = b(i) / m_A(i,i);
    }
    
    // Apply pivot matrix
    pivot(x);
}

//==============================================================================

template <typename Real>
void QRP<Real>::solveATA(Vector<Real> &x, Vector<Real>&b) const
{
    const int n = m_rank;
    
    // Apply the transpose of the pivot matrix
    pivot(b, true);

    // Forward-substitution with the implicit transpose of R (stored in A),
    // result stored in b
    for (int i = 0; i < n; ++i) {
        b(i) /= m_A(i,i);
        for (int j = i+1; j < n; ++j)
            b(j) -= m_A(i,j) * b(i);
    }
    
    // Back-subsitution
    for (int i = n-1; i >= 0; --i) {
        for (int j = i+1; j < n; ++j)
            b(i) -= m_A(i,j) * x(j);
        x(i) = b(i) / m_A(i,i);
    }
    
    // Apply pivot matrix
    pivot(x);
}

//==============================================================================

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_QRP_H


