#ifndef NUMERICS_LDLT_H
#define NUMERICS_LDLT_H

#include "Vector.h"
#include "Matrix.h"

#include <iostream>
#include <cstdlib>

namespace Mutation {
    namespace Numerics {


/**
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
     * Constructs a new LDLT decomposition using the new matrix.
     */
    template <typename E>
    void setMatrix(const SymMatExpr<T, E>& mat);
    
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
void LDLT<T>::setMatrix(const SymMatExpr<T, E>& mat)
{
    
    m_L = mat;
    const size_t n = m_L.rows();
    m_work.resize(n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j)
            m_work(j) = m_L(j,j) * m_L(i,j);
        
        for (int j = 0; j < i; ++j)
            m_L(i,i) -= m_work(j) * m_L(i,j);
        
        if (m_L(i,i) == static_cast<T>(0)) {
            std::cerr << "Calling LDLT decomposition with singular matrix!" << std::endl;
            exit(1);
        }
        
        for (int j = 0; j < i; ++j)
            for (int k = i+1; k < n; ++k)
                m_L(k,i) -= m_L(k,j) * m_work(j);
        
        for (int k = i+1; k < n; ++k)
            m_L(k,i) /= m_L(i,i);
    }
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
