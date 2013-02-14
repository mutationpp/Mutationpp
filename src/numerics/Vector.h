#ifndef NUMERICS_VECTOR_H
#define NUMERICS_VECTOR_H

#include <vector>
#include <cassert>
#include <cstring>
#include <iostream>

#include "Functors.h"
#include "ReferenceServer.h"

namespace Mutation {
    namespace Numerics {



template<typename T, typename E> class MatExpr;

/**
 * Represents an (index,value) pair from a vector expression.
 */
template <typename T>
class IVP : public std::pair<size_t, T>
{
public:
    using std::pair<size_t, T>::first;
    using std::pair<size_t, T>::second;
    
    IVP(const size_t& i, const T& v) 
        : std::pair<size_t, T>(i, v) { }
    
    const size_t& index() { return first; }
    const T& value() { return second; }
};

/**
 * Base class for all vector expressions.  Allows for static polymorphism 
 * between derived expression templates.
 */
template <typename T, typename Expr>
class VecExpr
{
public:
    /**
     * Returns derived class's size().
     */
    size_t size() const { return static_cast<const Expr&>(*this).size(); }

    /**
     * Returns derived class's operator() for writing.
     */
    T& operator()(const size_t i) {
        return static_cast<Expr&>(*this)(i); 
    }

    /**
     * Returns derived class's operator() for reading.
     */
    T operator()(const size_t i) const {
        return static_cast<const Expr&>(*this)(i); 
    }
    
    /**
     * Returns non-const reference to derived class.
     */
    operator Expr&() { return static_cast<Expr&>(*this); }
    
    /**
     * Returns const refernce to derived class.
     */
    operator const Expr&() const { return static_cast<const Expr&>(*this); }

    #define VEC_MAXMIN(__name__,__lg__)\
    IVP<T> __name__ () const {\
        const size_t n = static_cast<const Expr&>(*this).size();\
        T val = (n > 0 ? static_cast<const Expr&>(*this)(0) : T());\
        size_t k = (n > 0 ? 0 : -1);\
        for (size_t i = 1; i < n; ++i) {\
            if (static_cast<const Expr&>(*this)(i) __lg__ val) {\
                val = static_cast<const Expr&>(*this)(i);\
                k = i;\
            }\
        }\
        return IVP<T>(k, val);\
    }
    
    template <typename E>
    bool operator==(const VecExpr<T, E>& vec) {
        if (size() != vec.size())
            return false;
        for (int i = 0; i < size(); ++i)
            if (std::abs(static_cast<const Expr&>(*this)(i) - vec(i)) > 
                NumConst<T>::eps * 10.0)
                return false;
        return true;
    }
    
    template <typename E>
    bool operator!=(const VecExpr<T, E>& vec) {
        return !operator==(vec);
    }
    
    /**
     * Returns the pair (k, max), where Vector(k) = max is the maximum (ie: most 
     * positive) element in the vector.  If the vector size is zero, the pair
     * (-1, 0) is returned.
     */
    VEC_MAXMIN(max, >)
    
    /**
     * Returns the pair (k, min), where Vector(k) = min is the minimum (ie: most 
     * negative) element in the vector.  If the vector size is zero, the pair
     * (-1, 0) is returned.
     */
    VEC_MAXMIN(min, <)
    #undef VEC_MAXMIN
    
    /**
     * Returns the sum of the elements in the vector.
     */
    T sum() const { return accumulate(T(), Plus<T>()); }
    
    /**
     * Returns the 1-norm of the vector.
     */
    T norm1() const { return accumulate(T(), PlusAbs<T>()); }
    
    /**
     * Returns the sum of the squares of each element in the vector.
     */
    T sumOfSquares() const { return accumulate(T(), PlusSquare<T>()); }
    
    /**
     * Returns the 2-norm of the vector.
     */
    T norm2() const { return std::sqrt(sumOfSquares()); }
    
    /**
     * Returns the infinity norm of the vector.
     */
    T normInf() const { return accumulate(T(), MaxAbs<T>()); }
    
    /**
     * Returns the maximum (most positive) value in the vector.
     */
    T maxValue() const { return accumulate(T(), Max<T>()); }
    
    /**
     * Returns the minimum (most negative) value in the vector.
     */
    T minValue() const { return accumulate(T(), Min<T>()); }

protected:
    
    /**
     * Acts like the std::accumulate method.
     */
    template <typename BinOp>
    inline T accumulate(const T& init, const BinOp& op) const
    {
        T val = init;
        for (size_t i = 0; i < size(); ++i)
            val = op(val, static_cast<const Expr&>(*this)(i));
        return val;
    }
    
}; // class VecExpr


/**
 * Represents a slice out of a (possibly) larger Vector.  Created by accessing
 * a set of Vector elements with a start index and an end index.
 *
 * @see Vector<T>::operator(const int start, const int end)
 */
template <typename T, typename E>
class VectorSlice : public VecExpr<T, VectorSlice<T, E> >
{
public:

    VectorSlice() : mp_ref(NULL), m_start(0), m_size(0) { }    
    
    /**
     * Returns the size of this slice.
     */
    size_t size() const { return m_size; }
    
    /**
     * Accesses the i'th element for writing.  Bounds checking left to the
     * underlying Vector.
     */
    T& operator()(const size_t i) { 
        return static_cast<E&>(*mp_ref)(m_start + i); }
    
    /**
     * Accesses the i'th element for reading.  Bounds checking left to the
     * underlying Vector.
     */
    T operator()(const size_t i) const {
        return static_cast<const E&>(*mp_ref)(m_start + i);
    }
    
    #define VEC_SLICE_EQ_VEC(__op__,__assertion__)\
    template <typename T2, typename E2>\
    VectorSlice<T, E>& operator __op__ (const VecExpr<T2, E2>& vec) {\
        assert( size() == vec.size() );\
        const E2& v = vec;\
        for (size_t i = 0; i < size(); ++i) {\
            assert( __assertion__ );\
            static_cast<E&>(*mp_ref)(m_start + i) __op__ static_cast<T2>(v(i));\
        }\
        return *this;\
    }
    
    /**
     * Copies the value of the Vector expression vec to the underlying subvector
     * represented by this VectorSlice.  Assertion error thrown for different
     * sizes.
     */
    VEC_SLICE_EQ_VEC(=, true)
    
    /**
     * Adds the Vector expression vec to the underlying subvector represented
     * by this VectorSlice.  Assertion error thrown for different sizes.
     */
    VEC_SLICE_EQ_VEC(+=, true)
    
    /**
     * Subracts the Vector expression vec to the underlying subvector 
     * represented by this VectorSlice.  Assertion error thrown for different 
     * sizes.
     */
    VEC_SLICE_EQ_VEC(-=, true)
    
        /**
     * Multiplies the Vector expression with the underlying subvector
     * (element-wise) and stores solution in the underlying subvector.  
     * Assertion error thrown for different sizes.
     */
    VEC_SLICE_EQ_VEC(*=, true)
    
    /**
     * Divides the Vector by the Vector expression (element-wise) and stores
     * solution in the Vector.  Assertion error for divide by zero or
     * expressions with different sizes.
     */
    VEC_SLICE_EQ_VEC(/=, static_cast<T>(v(i)) != T())
    #undef VEC_SLICE_EQ_VEC
    
    #define VEC_SLICE_EQ_MAT(__op__)\
    template <typename T2, typename E2>\
    VectorSlice<T, E>& operator __op__ (const MatExpr<T2, E2>& mat) {\
        assert( (size() == mat.rows() && mat.cols() == 1) ||\
                (size() == mat.cols() && mat.rows() == 1) );\
        if (mat.rows() == 1) {\
            for (size_t i = 0; i < size(); ++i)\
                static_cast<E&>(*mp_ref)(m_start + i) __op__ mat(0,i);\
        } else {\
            for (size_t i = 0; i < size(); ++i)\
                static_cast<E&>(*mp_ref)(m_start + i) __op__ mat(i,0);\
        }\
        return *this;\
    }
    
    /**
     * Copies a Matrix expression (with at least one dimension equal to 1) to 
     * the underlying subvector represented by this VectorSlice.  Assertion
     * error for 2-dimensional matrix expressions or if the long dimension size
     * is not equal to the size of this slice.
     */
    VEC_SLICE_EQ_MAT(=)
    
    /**
     * Adds a Matrix expression (with at least one dimension equal to 1) to 
     * the underlying subvector represented by this VectorSlice.  Assertion
     * error for 2-dimensional matrix expressions or if the long dimension size
     * is not equal to the size of this slice.
     */
    VEC_SLICE_EQ_MAT(+=)
    
    /**
     * Subtracts a Matrix expression (with at least one dimension equal to 1) to 
     * the underlying subvector represented by this VectorSlice.  Assertion
     * error for 2-dimensional matrix expressions or if the long dimension size
     * is not equal to the size of this slice.
     */
    VEC_SLICE_EQ_MAT(-=)
    #undef VEC_SLICE_EQ_MAT
    
    /**
     * Sets the underlying Vector expression type and bounds for this slice.
     * Assertion errors thrown for bad indices or NULL pointer to expression.
     */
    VectorSlice<T, E>& setSlice(
        E *const p, const size_t start, const size_t end) {
        assert( p != NULL );
        assert( end <= static_cast<const E&>(*p).size() );
        assert( start < end );
        mp_ref = p;
        m_start = start;
        m_size = end - start;
        return *this;
    }
    
private:

    E* mp_ref;
    size_t m_start;
    size_t m_size;
    
}; // class VectorSlice


/**
 * Represents a numerical vector which can be manipulated via standard vector
 * algebra operations.
 */
template <typename T>
class Vector : public VecExpr<T, Vector<T> >
{
public:

    #define ASSERT_IN_RANGE(__index__)\
        assert( __index__ < size() );
    
    /**
     * Constructs a Vector with size n and intializes all values according to
     * init.
     */
    explicit Vector(size_t n = 0, const T& init = T()) : m_data(n, init) { }
    
    /**
     * Constructs a Vector from a c-array and size.
     */
    template <typename T2>
    Vector(const T2 *const p_data, const size_t n) 
        : m_data(p_data, p_data + n) { }
    
    /**
     * Copy constructor.
     */
    template <typename E>
    Vector(const VecExpr<T, E>& vec) {
        operator=(vec);
    }
    
    /**
     * Returns the size of this Vector.
     */
    size_t size() const { return m_data.size(); }
    
    /**
     * Returns the Vector element at index i for writing.  Assertion error
     * thrown for index out of Vector bounds.
     */
    T& operator()(const size_t i) {
        ASSERT_IN_RANGE(i);
        return static_cast<T&>(m_data[i]); 
    }
    
    /**
     * Returns the Vector element at index i for reading.  Assertion error
     * thrown for index out of Vector bounds.
     */
    T operator()(const size_t i) const {
        ASSERT_IN_RANGE(i);
        return m_data[i]; 
    }
    
    /**
     * Returns a VectorSlice for writing corresponding to the elements in 
     * indices [start, end).
     */
    VectorSlice<T, Vector<T> >& operator()(
        const size_t start, const size_t end) {
        return Mutation::Utilities::ReferenceServer<
            VectorSlice<T, Vector<T> >, 5>::serve().setSlice(this, start, end);
    }
    
    /**
     * Returns a VectorSlice for reading corresponding to the elements in 
     * indices [start, end).
     */
    const VectorSlice<T, Vector<T> >& operator()(
        const size_t start, const size_t end) const {
        return Mutation::Utilities::ReferenceServer<
            VectorSlice<T, Vector<T> >, 5>::serve().setSlice(
                const_cast<Vector<T>*>(this), start, end);
    }
    
    #define VEC_EQ_SCALAR(__op__,__assertion__)\
    Vector<T>& operator __op__ (const T& scalar) {\
        assert( __assertion__ );\
        for (size_t i = 0; i < size(); ++i)\
            m_data[i] __op__ scalar;\
        return *this;\
    }
    
    /**
     * Sets all values of the vector equal to scalar.  Returns the updated
     * Vector.
     */
    VEC_EQ_SCALAR(=,true)
    
    /**
     * Adds scalar to each element of the vector. Returns the updated Vector.
     */
    VEC_EQ_SCALAR(+=,true)
    
    /**
     * Subtracts scalar from each element of the vector. Returns the updated 
     * Vector. 
     */
    VEC_EQ_SCALAR(-=,true)
    
    /**
     * Multiplies each element of the vector by scalar. Returns the updated 
     * Vector. 
     */
    VEC_EQ_SCALAR(*=,true)
    
    /**
     * Divides each element of the vector by a scalar. Returns the updated 
     * Vector.  Assertion error thrown for attempt to divide by zero.
     */
    VEC_EQ_SCALAR(/=,scalar != T())
    #undef VEC_EQ_SCALAR
    
    /**
     * Overloaded assignment operator.  Resizes the vector to match vec if 
     * necessary.
     */
    template <typename T2, typename E>
    Vector<T>& operator=(const VecExpr<T2, E>& vec) {
        const E& v = vec;
        m_data.resize(v.size());
        for (size_t i = 0; i < v.size(); ++i)
            m_data[i] = static_cast<T>(v(i));
        return *this;
    }
        
    #define VEC_EQ_VEC(__op__,__assertion__)\
    template <typename T2, typename E>\
    Vector<T>& operator __op__ (const VecExpr<T2, E>& vec) {\
        assert( vec.size() == size() );\
        const E& v = vec;\
        for (size_t i = 0; i < size(); ++i) {\
            assert( __assertion__ );\
            m_data[i] __op__ static_cast<T>(v(i));\
        }\
        return *this;\
    }
    
    /**
     * Adds the Vector expression to this Vector.  Assertion error thrown for 
     * Vectors with non-equal sizes.
     */
    VEC_EQ_VEC(+=, true)
    
    /**
     * Subtracts the Vector expression from this Vector. Assertion error thrown 
     * for Vectors with nonequal sizes.
     */
    VEC_EQ_VEC(-=, true)
    
    /**
     * Multiplies the Vector expression with this Vector (element-wise) and
     * stores solution in the Vector.
     */
    VEC_EQ_VEC(*=, true)
    
    /**
     * Divides the Vector by the Vector expression (element-wise) and stores
     * solution in the Vector.  Assertion error for divide by zero.
     */
    VEC_EQ_VEC(/=, static_cast<T>(v(i)) != T())
    #undef VEC_EQ_VEC
    
    /**
     * Normalizes the vector to have unit length.  The only exception here is
     * that if the vector equals zero, it remains a zero vector.
     */
    Vector<T>& normalize() {
        T norm = static_cast<VecExpr<T, Vector<T> >&>(*this).norm2();
        return (norm != T() ? operator/=(norm) : *this);
    }
    
    /**
     * Replaces this Vector with a vector of the natural log of each element and
     * returns the Vector.
     */
    Vector<T>& log() { return applyFunc(Log<T>()); }
    
    /**
     * Replaces this Vector with a vector of the exponential of each element and
     * returns the Vector.
     */
    Vector<T>& exp() { return applyFunc(Exp<T>()); }
    
    /**
     * Replaces this Vector with a vector of the square-roots of each element
     * and returns the Vector.
     */
    Vector<T>& sqrt() { return applyFunc(Sqrt<T>()); }
    
    /**
     * Swaps element i and j in the Vector.  Assertion error thrown for indices
     * out of Vector bounds [0, size()).
     */
    void swap(const size_t& i, const size_t& j) {   
        ASSERT_IN_RANGE(i);
        ASSERT_IN_RANGE(j);
        if (i != j) {
            T temp = m_data[i];
            m_data[i] = m_data[j];
            m_data[j] = temp;
        }
    }
    
    #undef ASSERT_IN_RANGE

private:

    /**
     * Replaces the vector with a vector in which a unary operator is applied to
     * all the elements of this vector.
     */
    template <typename UnaryOp>
    inline Vector<T>& applyFunc(const UnaryOp& op)
    {
        for (size_t i = 0; i < size(); ++i)
            m_data[i] = op(m_data[i]);
        return *this;
    }
    
private:

    std::vector<T> m_data;

}; // class Vector


// Binary operators for Vector expressions
#define VEC_BINARY_VEC(__name__, __op__,__assertion__)\
template <typename T, typename E1, typename E2>\
class __name__ : public VecExpr<T, __name__ <T, E1, E2> >\
{\
public:\
    __name__ (const VecExpr<T, E1>& u, const VecExpr<T, E2>& v)\
        : m_u(u), m_v(v)\
    {\
        assert( u.size() == v.size() );\
    }\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i) const {\
        assert( __assertion__ );\
        return m_u(i) __op__ m_v(i);\
    }\
private:\
    const E1& m_u;\
    const E2& m_v;\
};\
template <typename T, typename E1, typename E2>\
const __name__ <T, E1, E2> operator __op__ (\
    const VecExpr<T, E1>& u, const VecExpr<T, E2>& v)\
{\
    return __name__ <T, E1,E2>(u, v);\
}
VEC_BINARY_VEC(VecDiffVec,-,true)
VEC_BINARY_VEC(VecPlusVec,+,true)
VEC_BINARY_VEC(VecMultVec,*,true)
VEC_BINARY_VEC(VecDivVec,/,m_v(i) != T())
#undef VEC_BINARY_OP

// Binary operators between scalars and vectors (and vectors and scalars)
#define VEC_BINARY_SCALAR(__name__, __op__,__assertion__)\
template <typename T, typename E>\
class __name__ : public VecExpr<T, __name__ <T, E> >\
{\
public:\
    __name__(const T alpha, const VecExpr<T, E>& u)\
        : m_alpha(alpha), m_u(u)\
    { assert( __assertion__ );}\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i) const { return m_u(i) __op__ m_alpha; }\
private:\
    const T m_alpha;\
    const E& m_u;\
};\
template <typename T, typename E>\
const __name__ <T, E> operator __op__ (const VecExpr<T, E>& u, const T alpha) {\
    return __name__ <T, E>(alpha, u);\
}

#define SCALAR_BINARY_VEC(__name__, __op__)\
template <typename T, typename E>\
class __name__ : public VecExpr<T, __name__ <T, E> >\
{\
public:\
    __name__(const T alpha, const VecExpr<T, E>& u)\
        : m_alpha(alpha), m_u(u) { }\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i) const { return m_alpha __op__ m_u(i); }\
private:\
    const T m_alpha;\
    const E& m_u;\
};\
template <typename T, typename E>\
const __name__ <T, E> operator __op__ (const T alpha, const VecExpr<T, E>& u) {\
    return __name__<T, E>(alpha, u);\
}

VEC_BINARY_SCALAR(VecDiffScalar,-,true)
VEC_BINARY_SCALAR(VecPlusScalar,+,true)
VEC_BINARY_SCALAR(VecProdScalar,*,true)
VEC_BINARY_SCALAR(VecDivScalar,/,alpha != T())

SCALAR_BINARY_VEC(ScalarDiffVec,-)
SCALAR_BINARY_VEC(ScalarPlusVec,+)
SCALAR_BINARY_VEC(ScalarProdVec,*)
SCALAR_BINARY_VEC(ScalarDivVec,/)

#undef VEC_BINARY_SCALAR
#undef SCALAR_BINARY_VEC

/**
 * Computes the inner product of two Vector expressions.  Throws assertion error
 * for expressions with different sizes.
 */
template <typename T, typename E1, typename E2>
T dot(const VecExpr<T, E1>& u, const VecExpr<T, E2>& v)
{
    assert( u.size() == v.size() );
    T sum = T();
    for (size_t i = 0; i < u.size(); ++i)
        sum += u(i) * v(i);
    return sum;
}

#define VEC_APPLY_FUNC(__name__,__class__,__op__)\
template <typename T, typename E>\
class __class__ : public VecExpr<T, __class__ <T, E> >\
{\
public:\
    __class__ (const VecExpr<T, E>& u) : m_u(u) { }\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i) const { return __op__ ; }\
private:\
    const E& m_u;\
};\
template <typename T, typename E>\
__class__ <T, E> __name__ (const VecExpr<T, E>& u) {\
    return __class__ <T, E>(u);\
}

VEC_APPLY_FUNC(square, VecSquared, m_u(i) * m_u(i))
VEC_APPLY_FUNC(log, VecLog, std::log(m_u(i)))
VEC_APPLY_FUNC(exp, VecExp, std::exp(m_u(i)))
#undef VEC_APPLY_FUNC

/*template <typename T>
class OnesVec : public VecExpr<T, OnesVec<T> >
{
public:
    OnesVec(const size_t size) : m_size(size) { }
    size_t size() const { return m_size; }
    T operator()(const size_t i) const { return static_cast<T>(1); }
private:
    const size_t m_size;
};*/

template <typename T>
class ConstVec : public VecExpr<T, ConstVec<T> >
{
public:
    ConstVec(const T& alpha, const size_t size) 
        : m_alpha(alpha), m_size(size) { }
    size_t size() const { return m_size; }
    T operator()(const size_t i) const { return m_alpha; }
private:
    const T& m_alpha;
    const size_t m_size;
};

template <typename T>
ConstVec<T> ones(const size_t size) {
    return ConstVec<T>(static_cast<T>(1), size);
}


template <typename T, typename E1, typename E2>
class MaxVec : public VecExpr<T, MaxVec<T, E1, E2> >
{
public:
    MaxVec(const VecExpr<T, E1>& u, const VecExpr<T, E2>& v)
        : m_u(u), m_v(v)
    {
        assert( u.size() == v.size() );
    }
    size_t size() const { return m_u.size(); }
    T operator()(const size_t i) const { return std::max(m_u(i), m_v(i)); }
private:
    const E1& m_u;
    const E2& m_v;
};

template <typename T, typename E1>
class MaxVecScalar : public VecExpr<T, MaxVecScalar<T, E1> >
{
public:
    MaxVecScalar(const VecExpr<T, E1>& u, const T& alpha)
        : m_u(u), m_alpha(alpha)
    { }
    size_t size() const { return m_u.size(); }
    T operator()(const size_t i) const { return std::max(m_u(i), m_alpha); }
private:
    const E1& m_u;
    const T& m_alpha;
};

template <typename T, typename E1, typename E2>
MaxVec<T, E1, E2> max(const VecExpr<T, E1>& u, const VecExpr<T, E2>& v) {
    return MaxVec<T, E1, E2>(u, v);
}
template <typename T, typename E>
MaxVecScalar<T, E> max(
    const T& alpha, const VecExpr<T, E>& u) {
    return MaxVecScalar<T, E>(u, alpha);
}
template <typename T, typename E>
MaxVecScalar<T, E > max(
    const VecExpr<T, E>& u, const T& alpha) {
    return MaxVecScalar<T, E>(u, alpha);
}

template <typename T>
class VectorWrapper : public VecExpr<T, VectorWrapper<T> >
{
public:
    VectorWrapper(T *const p_data, const size_t size) 
        : mp_data(p_data), m_size(size)
    {
        assert( p_data != NULL );
    }
    size_t size() const { return m_size; }
    T operator()(const size_t i) const { 
        assert( i < m_size );
        return mp_data[i]; 
    }
    T& operator()(const size_t i) {
        assert( i < m_size );
        return mp_data[i];
    }
    const VectorSlice<T, VectorWrapper<T> >& operator()(
        const size_t start, const size_t end) const {
        return Mutation::Utilities::ReferenceServer<
            VectorSlice<T, VectorWrapper<T> >, 5>::serve().setSlice(
                const_cast<VectorWrapper<T>*>(this), start, end);
    }
    
    template <typename E>
    VectorWrapper& operator=(const VecExpr<T, E>& vec) {
        assert( vec.size() == size() );
        for (int i = 0; i < size(); ++i)
            (*this)(i) = vec(i);
    }
private:
    T *const mp_data;
    const size_t m_size;
};

template <typename T>
VectorWrapper<T> asVector(T *const p_data, const size_t size) {
    return VectorWrapper<T>(p_data, size);
}

template <typename T>
const VectorWrapper<T> asVector(const T *const p_data, const size_t size) {
    return VectorWrapper<T>(const_cast<T *const>(p_data), size);
}

/**
 * Sets format flags for ostream depending on the data in array.
 */
template <typename T, typename ArrayType>
void setFancyFormat(std::ostream& os, const ArrayType& array, size_t& width, 
    std::ios::fmtflags& flags)
{
    const size_t max_width = 11;
    const size_t n = array.size();
    
    bool all_ints = true;
    int sig_upper = -1;
    int sig_lower =  0;
    int digit;
    T x, intpart;
    for (size_t i = 0; i < n; ++i) {
        x = std::abs(array(i));
        
        if (std::modf(x, &intpart) != static_cast<T>(0))
            all_ints = false;
        
        digit = (x != T() ? std::floor(std::log10(x)) : 0);
        sig_upper = std::max(sig_upper, digit);
        sig_lower = std::min(sig_lower, digit);
    }
    
    flags = os.flags();
    sig_upper += 1;    
    
    if (all_ints) {
        width = sig_upper + 2;
        if (width > max_width) {
            // Integers but too large to fit without using exponential form
            width = max_width;
            os.precision(max_width-8);
            os.setf(std::ios::scientific, std::ios::floatfield);
        } 
            // Otherwise, integers that fit in the maximum width
    } else {
        width = sig_upper - sig_lower + 5;
        if (width > max_width) {
            // Floating point numbers that must use exponential form to fit with
            // desired precision
            width = max_width;
            os.precision(max_width-8);
            os.setf(std::ios::scientific, std::ios::floatfield);
        } else {
            // Floating point numbers that can fit with desired precision
            os.precision(2-sig_lower);
            os.setf(std::ios::fixed, std::ios::floatfield);
        }
    }
    
    os.setf(std::ios::right, std::ios::adjustfield);
}

/**
 * Writes a vector expression to the output stream using fancy formating based
 * on the type of the vector.
 */
template <typename T, typename E> 
std::ostream& operator<<(std::ostream& os, const VecExpr<T, E>& v)
{
    const size_t n = v.size();
    
    std::ios::fmtflags flags;
    size_t width;
    
    setFancyFormat<T, Vector<T> >(os, v, width, flags);
    
    for (size_t i = 0; i < n; ++i)
        os << std::setw(width) << v(i) << std::endl;
    
    os.flags(flags);
    
    return os;
}

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_VECTOR_H
