#ifndef NUMERICS_PRECONDITIONER_H
#define NUMERICS_PRECONDITIONER_H

#include "Vector.h"
#include "Matrix.h"

namespace Numerics {

template <typename T>
class DiagonalPreconditioner
{
public:
    template <typename E>
    DiagonalPreconditioner(const MatExpr<T, E>& mat)
        : m_vec(mat.rows())
    {
        for (size_t i = 0; i < m_vec.size(); ++i)
            m_vec(i) = 1.0 / mat(i,i);
    }
    
    template <typename E>
    inline VecMultVec<T, Vector<T>, E> solve(VecExpr<T, E>& rhs) const {
        return rhs * m_vec;
    }
private:
    Vector<T> m_vec;   
};

template <typename T>
class IdentityPreconditioner
{
public:
    IdentityPreconditioner() { }
    
    template <typename E>
    inline VecExpr<T,E>& solve(VecExpr<T, E>& rhs) const {
        return rhs;
    }
};

} // namespace Numerics


#endif // NUMERICS_PRECONDITIONER_H
