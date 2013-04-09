/**
 * @file QR.h
 *
 * Implements the QR class which performs the QR decomposition.
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

//==============================================================================

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

//==============================================================================

template <typename Real>
QR<Real>::QR(const Matrix<Real> &A, const bool do_nothing)
    : m_A(A), m_betas(std::min(A.rows()-1,A.cols())), m_have_Q(false), 
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

//==============================================================================

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

//==============================================================================

template <typename Real>
const Matrix<Real> &QR<Real>::Q()
{
    assert( "Not implemented!" == 0 );
}

//==============================================================================

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

//==============================================================================

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

//==============================================================================

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

//==============================================================================

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_QR_H
