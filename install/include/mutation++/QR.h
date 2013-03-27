/**
 * @file QR.h
 *
 * Implements the QR and QRP classes which perform the QR decomposition without
 * and with column pivoting respectively. 
 *
 * @see class QR
 * @see class QRP
 *
 * @author J.B. Scoggins (jbscoggi@gmail.com)
 * @date   February 2, 2012
 */
 
#ifndef NUMERICS_QR_H
#define NUMERICS_QR_H

#include <cassert>
#include <vector>

#include "NumConst.h"
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {


/**
 * Computes the QR decomposition of a MxN matrix A where M >= N using 
 * Householder reflections.  The decomposition of A = QR is such that Q is MxM
 * orthogonal and R is MxN upper triangular.
 *
 * @todo Can move the house function into the code of QR and QRP directly and
 * act on the matrices themselves.  This would be more efficient because there
 * would no longer be the unnecessary copies back and forth.  Could probably
 * also use the norm vector computed in QRP and downdated with each transform
 * instead of recomputing the column norm for house...
 */
template <typename Real>
class QR
{
public:

    /**
     * Constructs the QR decomposition with column pivoting of A.
     */
    QR(const Matrix<Real>& A, bool do_nothing = false);
    
    /** 
     * Returns the NxN upper triangular matrix R.
     */
    const Matrix<Real> &R();
    
    /**
     * Returns the MxM orthogonal matrix Q. (TBD: just returns Q = 0 for now)
     *
     * @todo Make this function work if needed.
     */
    const Matrix<Real> &Q();
    
    /**
     * Computes Q'b.
     */
    void qTransposeB(Vector<Real> &b) const;
    
    /**
     * Solves Ax = b via Rx = Q'b.
     */
    void solve(Vector<Real> &x, Vector<Real> &b) const;
    
    /**
     * Solves the system A'Ax = b via R'Rx = b.
     */
    void solveATA(Vector<Real> &x, Vector<Real>&b) const;

protected:

    Matrix<Real> m_A;
    Matrix<Real> m_R;
    Matrix<Real> m_Q;
    
    bool m_have_Q;
    bool m_have_R;
    
    Vector<Real> m_betas;
};

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

/**
 * Computes the Householder vector associated with a given vector expression.
 * This algorithm is based on Algorithm 5.1.1 in "Matrix Computations" by Golub
 * and Van Loan.
 */
template<typename T, typename E>
void house(VecExpr<T, E>& v, T& beta)
{
    const size_t n = v.size();
    size_t i;
    
    // First compute ||v(2:n)||^2
    T sigma = NumConst<T>::zero;
    for (i = 1; i < n; ++i)
        sigma += v(i) * v(i);
    
    // Normalize v such that v(1) = 1    
    T v1 = v(0);
    if (sigma <= NumConst<T>::eps * abs(v1)) {
        v(0) = NumConst<T>::one;
        beta = NumConst<T>::zero;
    } else {
        T mu = sqrt(v1 * v1 + sigma);
        if (v1 <= NumConst<T>::zero)
            v1 -= mu;
        else
            v1 = -sigma / (v1 + mu);
        
        beta = NumConst<T>::two * v1 * v1 / (sigma + v1 * v1);
        
        v(0) = NumConst<T>::one;
        for (i = 1; i < n; ++i)
            v(i) /= v1;
    }
}

template <typename Real>
QR<Real>::QR(const Matrix<Real> &A, const bool do_nothing)
    : m_A(A), m_betas(min(A.rows()-1,A.cols())), m_have_Q(false), 
      m_have_R(false)
{
    // This modification allows the QRP class to inherit QR functionality,
    // without it, all of the common code between the two decompositions would
    // have to be repeated
    if (do_nothing)
        return;
    
    const size_t m = m_A.rows();
    const size_t n = m_A.cols();
    
    assert( m >= n );
    
    Vector<Real> u(m);
    Vector<Real> v(n);
    
    for (size_t i = 0; i < m_betas.size(); ++i) {
        // Compute Householder vector        
        house(u(i,m) = m_A(i,m,i,i+1), m_betas(i));
        
        // Update A
        v(i,n) = m_betas(i) * (u(i,m) * m_A(i,m,i,n));
        m_A(i,i) -= v(i);
        
        if (i+1 < n)
            m_A(i,m,i+1,n) -= rank1(u(i,m), v(i+1,n));
        
        // Store the Householder vector in A
        m_A(i+1,m,i,i+1) = u(i+1,m);
    }
}

template <typename Real>
const Matrix<Real> &QR<Real>::R()
{
    if (m_have_R)
        return m_R;
    
    m_R = Matrix<Real>(m_A.cols(), m_A.cols());
    
    for (int i = 0; i < m_R.cols(); ++i)
        for (int j = i; j < m_R.cols(); ++j)
            m_R(i,j) = m_A(i,j);
    
    m_have_R = true;
    return m_R;
}

template <typename Real>
const Matrix<Real> &QR<Real>::Q()
{
    assert( "Not implemented!" == 0 );
}

template <typename Real>
void QR<Real>::qTransposeB(Vector<Real> &b) const
{
    const size_t m = m_A.rows();
    const size_t n = m_A.cols();

    Vector<Real> u(m);
    for (size_t i = 0; i < m_betas.size(); ++i) {
        if (m_betas(i) > NumConst<Real>::zero) {
            u(i) = NumConst<Real>::one;
            u(i+1,m) = m_A(i+1,m,i,i+1);
            b(i,m) -= (m_betas(i) * dot(b(i,m), u(i,m))) * u(i,m);
        }
    }
}

template <typename Real>
void QR<Real>::solve(Vector<Real> &x, Vector<Real> &b) const
{
    // Update Q'b
    qTransposeB(b);
    
    // Back-substitution
    for (int i = m_A.cols()-1; i >= 0; --i) {
        for (int j = i+1; j < m_A.cols(); ++j)
            b(i) -= m_A(i,j) * x(j);
        x(i) = b(i) / m_A(i,i);
    }
}

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

template <typename Real>
void QR<Real>::solveATA(Vector<Real> &x, Vector<Real>&b) const
{
    const size_t n = m_A.cols();

    // Forward-substitution with the implicit transpose of R (stored in A),
    // result stored in b
    for (int i = 0; i < n; ++i) {
        b(i) /= m_A(i,i);
        for (int j = i+1; j < n; ++j)
            b(j) -= m_A(i,j) * b(i);
    }
    
    // Back-subsitution
    for (int i = n-1; i != 0; --i) {
        for (int j = i+1; j < n; ++j)
            b(i) -= m_A(i,j) * x(j);
        x(i) = b(i) / m_A(i,i);
    }
}

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

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_QR_H
