/**
 * @file Matrix.h
 *
 * @brief Defines the Matrix types and functionality.
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

#ifndef NUMERICS_MATRIX_H
#define NUMERICS_MATRIX_H

#include <iostream>
#include <iomanip>

#include "Vector.h"

namespace Mutation {
    namespace Numerics {


#define ASSERT_IN_RANGE(__low__,__index__,__high__)\
    assert( __index__ < __high__ );

/**
 * Base class for all Matrix expression types.  Allows static polymorphism
 * amongst the different expression types.
 */
template <typename T, typename E>
class MatExpr
{
public:

    /**
     * Returns derived class's rows().
     */
    size_t rows() const { return static_cast<const E&>(*this).rows(); }
    
    /**
     * Returns derived class's rows().
     */
    size_t cols() const { return static_cast<const E&>(*this).cols(); }
    
    /**
     * Returns derived class's size().
     */
    size_t size() const { return static_cast<const E&>(*this).size(); }
    
    /**
     * Returns derived matrix expression's operator(i).
     */
    T operator()(const size_t i) const {
        return static_cast<const E&>(*this)(i); 
    }
    
    /**
     * Returns derived matrix expression's operator(i).
     */
    T& operator()(const size_t i) {
        return static_cast<E&>(*this)(i); 
    }

    /**
     * Returns derived matrix expression's operator(i,j).
     */
    T operator()(const size_t i, const size_t j) const {
        return static_cast<const E&>(*this)(i, j); 
    }
    
    /**
     * Returns derived matrix expression's operator(i,j).
     */
    T& operator()(const size_t i, const size_t j) {
        return static_cast<E&>(*this)(i, j); 
    }
    
    /**
     * Returns non-const reference to derived class.
     */
    operator E&() { return static_cast<E&>(*this); }
    
    /**
     * Returns const refernce to derived class.
     */
    operator const E&() const { return static_cast<const E&>(*this); }
    
    /**
     * Returns the maximum element in this matrix.
     */
    T max() const {
        T maxval = (*this)(0);
        for (size_t i = 1; i < size(); ++i)
            maxval = std::max(maxval, (*this)(i));
        return maxval;
    }
    
}; // class MatExpr

/**
 * Base class for all Symmetric Matrix expressions.  Just used to enforce 
 * the correct operations between Symmetric and Non-Symmetric matrices.
 */
template <typename T, typename E>
class SymMatExpr : public MatExpr<T, E> { };


// Define some special matrices with delayed element evaluation (ie: a statement
// like Matrix<double> m = Identity<double>(10) will be evaluated in the Matrix
// constructor without creating a temporary matrix.
#define SPECIAL_MATRIX(__type__,__name__,__element__)\
template <typename T>\
class __name__ : public __type__ <T, __name__ <T> >\
{\
public:\
    __name__ (const size_t size)\
        : m_rows(size), m_cols(size) { }\
    __name__ (const size_t rows, const size_t cols)\
        : m_rows(rows), m_cols(cols) { }\
    size_t rows() const { return m_rows; }\
    size_t cols() const { return m_cols; }\
    T operator()(const size_t i, const size_t j) const {\
        return __element__ ;\
    }\
private:\
    const size_t m_rows;\
    const size_t m_cols;\
};

/**
 * Identity matrix.
 */
SPECIAL_MATRIX(SymMatExpr, Identity, static_cast<T>(i == j ? 1 : 0))

/**
 * Lehmer matrix.
 */
SPECIAL_MATRIX(SymMatExpr,
    Lehmer, static_cast<T>(std::min(i,j)+1) / static_cast<T>(std::max(i,j)+1))

/**
 * minij matrix.
 */
SPECIAL_MATRIX(SymMatExpr, MinIJ, static_cast<T>(std::min(i+1,j+1)))

/**
 * Zero matrix.
 */
SPECIAL_MATRIX(SymMatExpr, ZeroMat, static_cast<T>(0))
#undef SPECIAL_MATRIX

template <typename T>
ZeroMat<T> zeros(const size_t rows, const size_t cols) {
    return ZeroMat<T>(rows, cols);
}

/**
 * Represents a slice out of a (possibly) larger Matrix.  Created by accessing
 * a set of Matrix elements with a start index and an end index.
 *
 * @see Matrix<T>::operator(const size_t istart, const size_t iend, 
 *                          const size_t jstart, const size_t jend)
 */
template <typename T, typename E>
class MatrixSlice : public MatExpr<T, MatrixSlice<T, E> >
{
public:

    MatrixSlice() 
        : mp_ref(NULL), m_istart(0), m_jstart(0), m_rows(0), m_cols(0) { }

    /**
     * Returns number of rows of this slice.
     */
    size_t rows() const { return m_rows; }
    
    /**
     * Returns number of columns of this slice.
     */
    size_t cols() const { return m_cols; }
    
    /**
     * Returns the size of this slice.
     */
    size_t size() const { return m_rows * m_cols; }
    
    /**
     * Accesses the i'th element for writing.  Bounds checking left to the
     * underlying Vector.
     */
    T& operator()(const size_t i) { return (*this)(i / m_cols, i % m_cols); }
    
    /**
     * Accesses the i'th element for reading.  Bounds checking left to the
     * underlying Vector.
     */
    T operator()(const size_t i) const { 
        return (*this)(i / m_cols, i % m_cols);
    }
    
    /**
     * Returns the (i,j) element for writing.
     */
    T& operator()(const size_t i, const size_t j) {
        return static_cast<E&>(*mp_ref)(m_istart + i, m_jstart + j); 
    }

    /**
     * Returns the (i,j) element for reading.
     */
    T operator()(const size_t i, const size_t j) const {
        return static_cast<const E&>(*mp_ref)(m_istart + i, m_jstart + j); 
    }
    
    #define MAT_SLICE_EQ_SCALAR(__op__,__assertion__)\
    MatrixSlice<T, E>& operator __op__ (const T& scalar) {\
        assert( __assertion__ );\
        for (size_t i = 0; i < m_rows; ++i)\
            for (size_t j = 0; j < m_cols; ++j)\
                static_cast<E&>(*mp_ref)(m_istart + i, m_jstart + j) __op__ \
                    scalar;\
        return *this;\
    }
    
    MAT_SLICE_EQ_SCALAR(=,true)
    MAT_SLICE_EQ_SCALAR(+=,true)
    MAT_SLICE_EQ_SCALAR(-=,true)
    MAT_SLICE_EQ_SCALAR(*=,true)
    MAT_SLICE_EQ_SCALAR(/=,scalar != T())
    #undef MAT_SLICE_EQ_SCALAR
    
    /**
     * Copies the value of the Matrix expression mat to the underlying submatrix
     * represented by this MatrixSlice.  Assertion error thrown for different
     * sizes.
     */
    #define MAT_SLICE_EQ_MAT(__op__)\
    template <typename T2, typename E2>\
    MatrixSlice<T, E>& operator __op__ (const MatExpr<T2, E2>& mat) {\
        assert( m_rows == mat.rows() );\
        assert( m_cols == mat.cols() );\
        for (size_t i = 0; i < m_rows; ++i)\
            for (size_t j = 0; j < m_cols; ++j)\
                static_cast<E&>(*mp_ref)(m_istart + i, m_jstart + j) __op__ \
                    static_cast<T>(static_cast<const E2&>(mat)(i,j));\
        return *this;\
    }
    
    MatrixSlice<T, E>& operator=(const MatrixSlice<T, E>& mat) {
        assert( m_rows == mat.rows() );
        assert( m_cols == mat.cols() );
        const E& m = mat;
        for (size_t i = 0; i < m_rows; ++i)
            for (size_t j = 0; j < m_cols; ++j)
                static_cast<E&>(*mp_ref)(m_istart + i, m_jstart + j) =
                    static_cast<T>(m(i,j));
        return *this;
    }
    
    MAT_SLICE_EQ_MAT(=)
    MAT_SLICE_EQ_MAT(+=)
    MAT_SLICE_EQ_MAT(-=)
    #undef MAT_SLICE_EQ_MAT
    
    #define MAT_SLICE_EQ_VEC(__op__)\
    template <typename T2, typename E2>\
    MatrixSlice<T, E>& operator __op__ (const VecExpr<T2, E2>& vec) {\
        assert( (m_rows == vec.size() && m_cols == 1) ||\
                (m_cols == vec.size() && m_rows == 1) );\
        if (m_rows == 1) {\
            for (size_t i = 0; i < m_cols; ++i)\
                static_cast<E&>(*mp_ref)(m_istart, m_jstart + i) __op__ \
                    static_cast<T>(vec(i));\
        } else {\
            for (size_t i = 0; i < m_rows; ++i)\
                static_cast<E&>(*mp_ref)(m_istart + i, m_jstart) __op__ \
                    static_cast<T>(vec(i));\
        }\
        return *this;\
    }
    
    MAT_SLICE_EQ_VEC(=)
    MAT_SLICE_EQ_VEC(+=)
    MAT_SLICE_EQ_VEC(-=)
    #undef MAT_SLICE_EQ_VEC
    
    MatrixSlice<T, E>& setSlice(
        E *const p_ref, const size_t istart, const size_t iend, 
        const size_t jstart, const size_t jend) {
        assert( 0 <= istart );
        assert( 0 <= jstart );
        assert( istart < iend );
        assert( jstart < jend );
        assert( iend <= static_cast<E&>(*p_ref).rows() );
        assert( jend <= static_cast<E&>(*p_ref).cols() );
        mp_ref   = p_ref;
        m_istart = istart;
        m_jstart = jstart;
        m_rows   = iend - istart;
        m_cols   = jend - jstart;
        return *this;
    }
    
private:

    E* mp_ref;
    size_t m_istart;
    size_t m_jstart;
    size_t m_rows;
    size_t m_cols;
    
}; // class VectorSlice/**

/*
 * Represents a numeric matrix which can be manipulated by standard mathematical
 * operators for matrices.
 */
template <typename T>
class Matrix : public MatExpr<T, Matrix<T> >
{
public:
        
    /**
     * Empty constructor.
     */
    Matrix() : m_rows(0), m_cols(0), m_data(0) { }

    /**
     * Constructs an m by n matrix with all elements equal to init.
     */
    Matrix(const size_t m, const size_t n, const T& init = T())
        : m_rows(m), m_cols(n), m_data(m * n, init) { }
    
    /**
     * Copy constructor.
     */
    template <typename T2, typename E>
    Matrix(const MatExpr<T2, E>& mat) { operator=(mat); }
    
    /**
     * Returns number of rows of this matrix.
     */
    size_t rows() const { return m_rows; }
    
    /**
     * Returns number of columns of this matrix.
     */
    size_t cols() const { return m_cols; }
    
    /**
     * Returns the actual size of the underlying storage array.
     */
    size_t size() const { return m_data.size(); }
    
    /**
     * Returns the element indexed at i in the underlying storage array.
     */
    T& operator()(const size_t i) {
        ASSERT_IN_RANGE(0, i, m_data.size())
        return m_data[i];
    }
    
    /**
     * Returns the element indexed at i in the underlying storage array.
     */
    T operator()(const size_t i) const { 
        ASSERT_IN_RANGE(0, i, m_data.size())
        return m_data[i];
    }
    
    /**
     * Accesses an element at (i,j) for writing.  Assertion error thrown for 
     * indices out of bounds.
     */
    T& operator()(const size_t i, const size_t j) {
        ASSERT_IN_RANGE(0, i, m_rows);
        ASSERT_IN_RANGE(0, j, m_cols);
        return m_data[i * m_cols + j];
    }

    /**
     * Accesses an element at (i,j) for reading.  Assertion error thrown for 
     * indices out of bounds.
     */
    T operator()(const size_t i, const size_t j) const {
        ASSERT_IN_RANGE(0, i, m_rows);
        ASSERT_IN_RANGE(0, j, m_cols);
        return m_data[i * m_cols + j];
    }
    
    /**
     * Accesses a sub matrix for writing with rows in the range [i1,i2) and 
     * columns in the range [j1,j2).  
     */
    MatrixSlice<T, Matrix<T> >& operator()(
        const size_t i1, const size_t i2, const size_t j1, const size_t j2) {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i1, i2, j1, j2);
    }
    
    /**
     * Accesses a sub matrix for reading with rows in the range [i1,i2) and 
     * columns in the range [j1,j2).  
     */
    const MatrixSlice<T, Matrix<T> >& operator()(
        const size_t i1, const size_t i2, const size_t j1, const size_t j2)
        const {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i1, i2, j1, j2);
    }
    
    MatrixSlice<T, Matrix<T> >& row(const size_t i) {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i, i+1, 0, cols());
    }
    
    MatrixSlice<T, Matrix<T> >& col(const size_t i) {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), 0, rows(), i, i+1);
    }
    
    const MatrixSlice<T, Matrix<T> >& row(const size_t i) const {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i, i+1, 0, cols());
    }
    
    const MatrixSlice<T, Matrix<T> >& col(const size_t i) const {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), 0, rows(), i, i+1);
    }
    
    /**
     * Assign contents of a Matrix expression to this Matrix.  Resize if
     * necessary.
     */
    template <typename T2, typename E>
    Matrix<T>& operator=(const MatExpr<T2, E>& mat) {
        const E& m = mat;
        m_rows = m.rows();
        m_cols = m.cols();
        m_data.resize(m_rows * m_cols);
        for (size_t i = 0; i < m_rows; ++i)
            for (size_t j = 0; j < m_cols; ++j)
                (*this)(i,j) = static_cast<T>(m(i,j));
        return *this;
    }
    
    #define MAT_EQ_SCALAR(__op__,__assertion__)\
    Matrix<T>& operator __op__ (const T& scalar) {\
        assert( __assertion__ );\
        for (size_t i = 0; i < m_data.size(); ++i)\
            m_data[i] __op__ scalar;\
        return *this;\
    }
    
    /**
     * Sets all values of the matrix equal to scalar.  Returns the updated
     * Matrix.
     */
    MAT_EQ_SCALAR(=,true)
    
    /**
     * Adds scalar to each element of the matrix. Returns the updated Matrix.
     */
    MAT_EQ_SCALAR(+=,true)
    
    /**
     * Subtracts scalar from each element of the matrix. Returns the updated 
     * Matrix. 
     */
    MAT_EQ_SCALAR(-=,true)
    
    /**
     * Multiplies each element of the matrix by scalar. Returns the updated 
     * Matrix. 
     */
    MAT_EQ_SCALAR(*=,true)
    
    /**
     * Divides each element of the matrix by a scalar. Returns the updated 
     * Matrix.  Assertion error thrown for attempt to divide by zero.
     */
    MAT_EQ_SCALAR(/=,scalar != T())
    #undef MAT_EQ_SCALAR
    
    #define MAT_EQ_MAT(__op__)\
    template <typename T2, typename E>\
    Matrix<T>& operator __op__ (const MatExpr<T2, E>& mat) {\
        assert( mat.rows() == rows() );\
        assert( mat.cols() == cols() );\
        const E& m = mat;\
        for (size_t i = 0; i < rows(); ++i)\
            for (size_t j = 0; j < cols(); ++j)\
                (*this)(i,j) __op__ static_cast<T>(m(i,j));\
        return *this;\
    }
    
    /**
     * Adds the Matrix expression to this Matrix.  Assertion error thrown for 
     * Matrices with non-equal sizes.
     */
    MAT_EQ_MAT(+=)
    
    /**
     * Subtracts the Matrix expression from this Matrix. Assertion error thrown 
     * for Matrices with nonequal sizes.
     */
    MAT_EQ_MAT(-=)
    #undef MAT_EQ_MAT
    
    /**
     * Returns the transpose of the matrix.  Note that this matrix is not 
     * affected, instead a copy is created and transposed.
     */
    inline Matrix<T> transpose() const {
        Matrix<T> mat(cols(), rows());
        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                mat(j,i) = (*this)(i,j);
        return mat;
    }
    
    /**
     * Swaps the rows i and j.
     */
    inline void swapRows(const size_t i, const size_t j) {
        ASSERT_IN_RANGE(0, i, rows())
        ASSERT_IN_RANGE(0, j, rows())
        T temp;
        for (size_t k = 0; k < cols(); ++k) {
            temp = (*this)(i,k);
            (*this)(i,k) = (*this)(j,k);
            (*this)(j,k) = temp;
        }
    }
    
    /**
     * Swaps the columns i and j.
     */
    inline void swapCols(const size_t i, const size_t j) {
        ASSERT_IN_RANGE(0, i, cols())
        ASSERT_IN_RANGE(0, j, cols())
        T temp;
        for (size_t k = 0; k < rows(); ++k) {
            temp = (*this)(k,i);
            (*this)(k,i) = (*this)(k,j);
            (*this)(k,j) = temp;
        }
    }

private:

    typedef Mutation::Utilities::ReferenceServer<
        MatrixSlice<T, Matrix<T> >, 5> SliceServer;

    size_t m_rows;
    size_t m_cols;
    std::vector<T> m_data;
    
}; // class Matrix

/**
 * Represents a symmetric matrix which can be manipulated by standard
 * mathematical operators for symmetric matrices.
 */
template <typename T>
class SymmetricMatrix : public SymMatExpr<T, SymmetricMatrix<T> >
{
public:
       
    /**
     * Empty constructor.
     */
    SymmetricMatrix() : m_cols(0), m_data(0) { }

    /**
     * Constructs an m by n matrix with all elements equal to init.
     */
    SymmetricMatrix(const size_t n, const T& init = T())
        : m_cols(n), m_data((n*n - n) / 2 + n, init)
    {
        updateDiagonalIndices();
    }
    
    /**
     * Copy constructor.
     */
    template <typename T2, typename E>
    SymmetricMatrix(const SymMatExpr<T2, E>& mat) { operator=(mat); }
    
    /**
     * Returns number of rows of this matrix.
     */
    size_t rows() const { return m_cols; }
    
    /**
     * Returns number of columns of this matrix.
     */
    size_t cols() const { return m_cols; }
    
    /**
     * Returns the actual size of the underlying storage array.
     */
    size_t size() const { return m_data.size(); }
    
    /**
     * Returns the element indexed at i in the underlying storage array.
     */
    const T& operator()(const size_t i) const {
        ASSERT_IN_RANGE(0, i, m_data.size())
        return m_data[i];
    }
    
    T& operator()(const size_t i) { 
        ASSERT_IN_RANGE(0, i, m_data.size())
        return m_data[i];
    }
    
    /**
     * Accesses an element at (i,j) for writing.  Assertion error thrown for 
     * indices out of bounds.
     */
    T& operator()(const size_t i, const size_t j) {
        ASSERT_IN_RANGE(0, i, m_cols);
        ASSERT_IN_RANGE(0, j, m_cols);
        //return m_data[( i <= j ? i*m_cols - (i*i+i)/2 + j :
        //                         j*m_cols - (j*j+j)/2 + i )];
        return m_data[i <= j ? m_ii[i]+(j-i) : m_ii[j]+(i-j)];
    }

    /**
     * Accesses an element at (i,j) for reading.  Assertion error thrown for 
     * indices out of bounds.
     */
    T operator()(const size_t i, const size_t j) const {
        ASSERT_IN_RANGE(0, i, m_cols);
        ASSERT_IN_RANGE(0, j, m_cols);
        //return m_data[( i <= j ? i*m_cols - (i*i+i)/2 + j :
        //                         j*m_cols - (j*j+j)/2 + i )];
        return m_data[i <= j ? m_ii[i]+(j-i) : m_ii[j]+(i-j)];
    }
    
    /**
     * Accesses a sub matrix for writing with rows in the range [i1,i2) and 
     * columns in the range [j1,j2).  
     */
    MatrixSlice<T, Matrix<T> >& operator()(
        const size_t i1, const size_t i2, const size_t j1, const size_t j2) {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i1, i2, j1, j2);
    }
    
    /**
     * Accesses a sub matrix for reading with rows in the range [i1,i2) and 
     * columns in the range [j1,j2).  
     */
    const MatrixSlice<T, Matrix<T> >& operator()(
        const size_t i1, const size_t i2, const size_t j1, const size_t j2)
        const {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i1, i2, j1, j2);
    }
    
    MatrixSlice<T, SymmetricMatrix<T> >& row(const size_t i) {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), i, i+1, 0, cols());
    }
    
    MatrixSlice<T, SymmetricMatrix<T> >& col(const size_t i) {
        return SliceServer::serve().setSlice(
            const_cast<Matrix<T>*>(this), 0, rows(), i, i+1);
    }
    
    const MatrixSlice<T, SymmetricMatrix<T> >& row(const size_t i) const {
        return SliceServer::serve().setSlice(
            const_cast<SymmetricMatrix<T>*>(this), i, i+1, 0, cols());
    }
    
    const MatrixSlice<T, SymmetricMatrix<T> >& col(const size_t i) const {
        return SliceServer::serve().setSlice(
            const_cast<SymmetricMatrix<T>*>(this), 0, rows(), i, i+1);
    }
    
    /**
     * Assign contents of a SymmetricMatrix expression to this Matrix.  Resize 
     * if necessary.
     */
    template <typename T2, typename E>
    SymmetricMatrix<T>& operator=(const SymMatExpr<T2, E>& mat) {
        const E& m = mat;
        m_cols = m.cols();
        m_data.resize((m_cols*m_cols - m_cols)/2 + m_cols);
        for (int i = 0; i < m_data.size(); ++i)
            m_data[i] = m(i);
        updateDiagonalIndices();
        return *this;
    }
    
    SymmetricMatrix<T>& operator=(const SymmetricMatrix<T>& mat) {
        m_cols = mat.cols();
        m_data = mat.m_data;
        m_ii   = mat.m_ii;
        return *this;
    }

    #define MAT_EQ_SCALAR(__op__,__assertion__)\
    SymmetricMatrix<T>& operator __op__ (const T& scalar) {\
        assert( __assertion__ );\
        for (size_t i = 0; i < m_data.size(); ++i)\
            m_data[i] __op__ scalar;\
        return *this;\
    }
    
    /**
     * Sets all values of the matrix equal to scalar.  Returns the updated
     * Matrix.
     */
    MAT_EQ_SCALAR(=,true)
    
    /**
     * Adds scalar to each element of the matrix. Returns the updated Matrix.
     */
    MAT_EQ_SCALAR(+=,true)
    
    /**
     * Subtracts scalar from each element of the matrix. Returns the updated 
     * Matrix. 
     */
    MAT_EQ_SCALAR(-=,true)
    
    /**
     * Multiplies each element of the matrix by scalar. Returns the updated 
     * Matrix. 
     */
    MAT_EQ_SCALAR(*=,true)
    
    /**
     * Divides each element of the matrix by a scalar. Returns the updated 
     * Matrix.  Assertion error thrown for attempt to divide by zero.
     */
    MAT_EQ_SCALAR(/=,scalar != T())
    #undef MAT_EQ_SCALAR
    
    #define MAT_EQ_MAT(__op__)\
    template <typename T2, typename E>\
    SymmetricMatrix<T>& operator __op__ (const SymMatExpr<T2, E>& mat) {\
        assert( mat.cols() == cols() );\
        const E& m = mat;\
        for (size_t i = 0; i < rows(); ++i)\
            for (size_t j = i; j < cols(); ++j)\
                (*this)(i,j) __op__ static_cast<T>(m(i,j));\
        return *this;\
    }
    
    /**
     * Adds the Matrix expression to this Matrix.  Assertion error thrown for 
     * Matrices with non-equal sizes.
     */
    MAT_EQ_MAT(+=)
    
    /**
     * Subtracts the Matrix expression from this Matrix. Assertion error thrown 
     * for Matrices with nonequal sizes.
     */
    MAT_EQ_MAT(-=)
    #undef MAT_EQ_MAT
    
private:

    typedef Mutation::Utilities::ReferenceServer<
        MatrixSlice<T, Matrix<T> >, 5> SliceServer;

    // Store indices to diagonal elements
    void updateDiagonalIndices() {
        m_ii.resize(m_cols);
        for (int i = 0; i < m_cols; ++i)
            m_ii[i] = i*(m_cols+1) - (i*i+i)/2;
    }

    size_t m_cols;
    std::vector<T> m_data;
    std::vector<size_t> m_ii;
    
}; // class SymmetricMatrix

template <typename T>
class SymMatWrapper : public SymMatExpr<T, SymMatWrapper<T> >
{
public:
    SymMatWrapper(T *const p_data, const size_t cols) 
        : m_cols(cols), mp_data(p_data)
    { }
    size_t rows() const { return m_cols; }
    size_t cols() const { return m_cols; }
    size_t size() const { return (m_cols*m_cols + m_cols)/2; }
    
    T& operator()(const size_t i) {
        ASSERT_IN_RANGE(0,i,size());
        return mp_data[i];
    }
    
    T operator()(const size_t i) const {
        ASSERT_IN_RANGE(0,i,size());
        return mp_data[i];
    }
    
    /**
     * Accesses an element at (i,j) for writing.  Assertion error thrown for 
     * indices out of bounds.
     */
    T& operator()(const size_t i, const size_t j) {
        ASSERT_IN_RANGE(0, i, m_cols);
        ASSERT_IN_RANGE(0, j, m_cols);
        return mp_data[( i <= j ? i*m_cols - (i*i+i)/2 + j :
                                  j*m_cols - (j*j+j)/2 + i )];
    }

    /**
     * Accesses an element at (i,j) for reading.  Assertion error thrown for 
     * indices out of bounds.
     */
    T operator()(const size_t i, const size_t j) const {
        ASSERT_IN_RANGE(0, i, m_cols);
        ASSERT_IN_RANGE(0, j, m_cols);
        return mp_data[( i <= j ? i*m_cols - (i*i+i)/2 + j :
                                  j*m_cols - (j*j+j)/2 + i )];
    }
private:
    size_t m_cols;
    T* mp_data;
};


template <typename T, typename E>
class NegateMat : public MatExpr<T, NegateMat<T, E> >
{
public:
    NegateMat(const MatExpr<T, E>& mat) : m_A(mat) { }
    size_t rows() const { return m_A.rows(); }
    size_t cols() const { return m_A.cols(); }
    T operator()(const size_t i, const size_t j) const {
        return -m_A(i,j);
    }
private:
    const E& m_A;
};

template <typename T, typename E>
NegateMat<T, E> operator-(const MatExpr<T, E>& mat) {
    return NegateMat<T, E>(mat);
}

// Binary operators for Matrix expressions
#define MAT_BINARY_MAT(__name__, __op__)\
template <typename T, typename E1, typename E2>\
class __name__ : public MatExpr<T, __name__ <T, E1, E2> >\
{\
public:\
    __name__ (const MatExpr<T, E1>& u, const MatExpr<T, E2>& v)\
        : m_u(u), m_v(v)\
    {\
        assert( u.rows() == v.rows() );\
        assert( u.cols() == v.cols() );\
    }\
    size_t rows() const { return m_u.rows(); }\
    size_t cols() const { return m_u.cols(); }\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i, const size_t j) const {\
        return m_u(i,j) __op__ m_v(i,j);\
    }\
private:\
    const E1& m_u;\
    const E2& m_v;\
};\
template <typename T, typename E1, typename E2>\
const __name__ <T, E1, E2> operator __op__ (\
    const MatExpr<T, E1>& u, const MatExpr<T, E2>& v)\
{\
    return __name__ <T, E1,E2>(u, v);\
}
MAT_BINARY_MAT(MatDiffMat,-)
MAT_BINARY_MAT(MatPlusMat,+)
#undef VEC_BINARY_OP

// Binary operators between scalars and matrices (and matrices and scalars)
#define MAT_BINARY_SCALAR(__name__, __op__,__assertion__)\
template <typename T, typename E>\
class __name__ : public MatExpr<T, __name__ <T, E> >\
{\
public:\
    __name__(const T alpha, const MatExpr<T, E>& u)\
        : m_alpha(alpha), m_u(u)\
    { assert( __assertion__ );}\
    size_t rows() const { return m_u.rows(); }\
    size_t cols() const { return m_u.cols(); }\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i, const size_t j) const {\
        return m_u(i,j) __op__ m_alpha;\
    }\
private:\
    const T m_alpha;\
    const E& m_u;\
};\
template <typename T, typename E>\
const __name__ <T, E> operator __op__ (const MatExpr<T, E>& u, const T alpha) {\
    return __name__ <T, E>(alpha, u);\
}

#define SCALAR_BINARY_MAT(__name__,__op__)\
template <typename T, typename E>\
class __name__ : public MatExpr<T, __name__ <T, E> >\
{\
public:\
    __name__(const T alpha, const MatExpr<T, E>& u)\
        : m_alpha(alpha), m_u(u) { }\
    size_t rows() const { return m_u.rows(); }\
    size_t cols() const { return m_u.cols(); }\
    size_t size() const { return m_u.size(); }\
    T operator()(const size_t i, const size_t j) const {\
        return m_alpha __op__ m_u(i,j);\
    }\
private:\
    const T m_alpha;\
    const E& m_u;\
};\
template <typename T, typename E>\
const __name__ <T, E> operator __op__ (const T alpha, const MatExpr<T, E>& u) {\
    return __name__<T, E>(alpha, u);\
}

MAT_BINARY_SCALAR(MatDiffScalar,-,true)
MAT_BINARY_SCALAR(MatPlusScalar,+,true)
MAT_BINARY_SCALAR(MatProdScalar,*,true)
MAT_BINARY_SCALAR(MatDivScalar,/,alpha != T())

SCALAR_BINARY_MAT(ScalarDiffMat,-)
SCALAR_BINARY_MAT(ScalarPlusMat,+)
SCALAR_BINARY_MAT(ScalarProdMat,*)

#undef VEC_BINARY_SCALAR
#undef SCALAR_BINARY_MAT

// Matrix times vector (returns VecExpr)
template <typename T, typename VE, typename ME>
class MatProdVec : public VecExpr<T, MatProdVec<T, VE, ME> >
{ 
public: 
    MatProdVec(const VecExpr<T, VE>& x, const MatExpr<T, ME>& A)
        : m_x(x), m_A(A)
    {
        assert( x.size() == A.cols() );
    }
    
    size_t size() const { return m_A.rows(); }
    
    T operator()(const size_t i) const {
        T sum = T();
        for (size_t k = 0; k < m_x.size(); ++k)
            sum += m_A(i,k) * m_x(k);
        return sum;
    }
private:
    const VE& m_x;
    const ME& m_A;
};

template <typename T, typename VE, typename ME>
MatProdVec<T, VE, ME> operator*(
    const MatExpr<T, ME>& A, const VecExpr<T, VE>& x) {
    return MatProdVec<T, VE, ME>(x, A);
}

// Vector times matrix (returns VecExpr) (potentially slow for large matrices)
template <typename T, typename VE, typename ME>
class VecProdMat : public VecExpr<T, VecProdMat<T, VE, ME> >
{ 
public: 
    VecProdMat(const VecExpr<T, VE>& x, const MatExpr<T, ME>& A)
        : m_x(x), m_A(A)
    {
        assert( x.size() == A.rows() );
    }
    
    size_t size() const { return m_A.cols(); }
    
    T operator()(const size_t i) const {
        T sum = T();
        for (size_t k = 0; k < m_x.size(); ++k)
            sum += m_A(k,i) * m_x(k);
        return sum;
    }
private:
    const VE& m_x;
    const ME& m_A;
};

template <typename T, typename VE, typename ME>
VecProdMat<T, VE, ME> operator*(
    const VecExpr<T, VE>& x, const MatExpr<T, ME>& A) {
    return VecProdMat<T, VE, ME>(x, A);
}

// Matrix times matrix (returns MatExpr)
template <typename T, typename E1, typename E2>
class MatProdMat : public MatExpr<T, MatProdMat<T, E1, E2> >
{
public: 
    MatProdMat(const MatExpr<T, E1>& A, const MatExpr<T, E2>& B)
        : m_A(A), m_B(B)
    {
        assert( m_A.cols() == m_B.rows() );
    }
    
    size_t rows() const { return m_A.rows(); }
    size_t cols() const { return m_B.cols(); }
    size_t size() const { return m_A.rows() * m_B.cols(); }
    
    /*T operator()(const size_t i) const {
        return (*this)(i / m_B.cols(), i % m_B.cols());
    }*/
    
    T operator()(const size_t i, const size_t j) const {
        T sum = T();
        for (size_t k = 0; k < m_A.cols(); ++k)
            sum += m_A(i,k) * m_B(k,j);
        return sum;
    }
private:
    const E1& m_A;
    const E2& m_B;
};

// Via matrix expression (ie delayed)
template <typename T, typename E1, typename E2>
MatProdMat<T, E1, E2> operator*(
    const MatExpr<T, E1>& A, const MatExpr<T, E2>& B) {
    return MatProdMat<T, E1, E2>(A, B);
}

// Constructs a Rank1 update matrix (outer product of two vectors)
template <typename T, typename V1, typename V2>
class Rank1 : public MatExpr<T, Rank1<T, V1, V2> >
{ 
public: 
    Rank1(const VecExpr<T, V1>& u, const VecExpr<T, V2>& v)
        : m_u(u), m_v(v)
    { }
    
    size_t rows() const { return m_u.size(); }
    size_t cols() const { return m_v.size(); }
    
    T operator()(const size_t i, const size_t j) const {
        return m_u(i) * m_v(j);
    }
private:
    const V1& m_u;
    const V2& m_v;
};

template <typename T, typename V1, typename V2>
Rank1<T, V1, V2> rank1(const VecExpr<T, V1>& u, const VecExpr<T, V2>& v) {
    return Rank1<T, V1, V2>(u, v);
}

/**
 * Writes the matrix expression to the output stream using fancy formatting
 * based on the type of the matrix.
 */
template <typename T, typename E>
std::ostream& operator<<(std::ostream& os, const MatExpr<T, E>& mat)
{
    const size_t nr = mat.rows();
    const size_t nc = mat.cols();
    
    std::ios::fmtflags flags;
    size_t width;
    
    setFancyFormat<T, Matrix<T> >(os, mat, width, flags);
    
    for (size_t i = 0; i < nr; ++i) {
        for (size_t j = 0; j < nc; ++j)
            os << std::setw(width) << mat(i,j);
        os << std::endl;
    }
    
    os.flags(flags);
    
    return os;
}

#undef ASSERT_IN_RANGE

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_MATRIX_H
