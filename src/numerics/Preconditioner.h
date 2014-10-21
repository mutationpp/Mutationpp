#ifndef NUMERICS_PRECONDITIONER_H
#define NUMERICS_PRECONDITIONER_H

#include "Vector.h"
#include "Matrix.h"

namespace Mutation {
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
            m_vec(i) = (mat(i,i) != 0.0 ? 1.0 / mat(i,i) : 1.0);
    }
    
    template <typename E>
    inline const VecMultVec<T, Vector<T>, E> solve(const VecExpr<T, E>& rhs) const {
        return m_vec * rhs;
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
    inline const VecExpr<T,E>& solve(const VecExpr<T, E>& rhs) const {
        return rhs;
    }
};

    } // namespace Numerics
} // namespace Mutation


#endif // NUMERICS_PRECONDITIONER_H
