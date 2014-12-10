/**
 * @file LeastSquares.h
 *
 * @brief Implements the class LeastSquares which solves overdetermined linear
 * systems.
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

#ifndef NUMERICS_LEASTSQUARES_H
#define NUMERICS_LEASTSQUARES_H

//#define LS_WARN_RANK_DEFICIENT

#include <cassert>
#include <iostream>

#include "QRP.h"
#include "SVD.h"

namespace Mutation {
    namespace Numerics {


/**
 * @ingroup dirsol
 * Computes the minimum-norm least-squares solution to an overdetermined system 
 * Ax = b by storing intermediate calculations so that multiple (x,b) pairs do 
 * not require refactoring the matrix A.
 */
template <typename Real>
class LeastSquares
{
public:
    
    /**
     * Factors the matrix A for preparation to solve a least-squares problem.
     * A is factored using the rank-revealing QR decomposition with column-
     * pivoting.  If it is determined that A is full rank, nothing else is done,
     * however if A is rank deficient, the SVD of the NxN R matrix is also 
     * computed and stored.
     */
    template <typename E>
    LeastSquares(const MatExpr<Real, E> &A)
        : m_rows(A.rows()), m_cols(A.cols()), mp_qrp(NULL), mp_svd(NULL)
    {
        assert( A.rows() >= A.cols() );
        
        mp_qrp = new QRP<Real>(A);
        
        m_full_rank = (mp_qrp->rank() == A.cols());
        
        if (!m_full_rank)
            mp_svd = new SVD<Real>(mp_qrp->R());
    }
    
    /**
     * Destructor.
     */
    ~LeastSquares()
    {
        delete mp_qrp; 
        mp_qrp = NULL;
        
        if (!m_full_rank) {
            delete mp_svd;
            mp_svd = NULL;
        }
    }
    
    /**
     * Returns the rank of the matrix A.
     */
    int rank() const {
        return mp_qrp->rank();
    }
    
    /**
     * Solves the overdetermined system Ax = b.  If A is full rank then the
     * solution is obtained via RP'x = Q'b where AP = QR.  If A is rank 
     * deficient, the solution is computed by USV'P'x = Q'b where R = USV' is
     * the singular value decomposition of R from the QR decomposition.
     */
    void solve(Vector<Real> &x, Vector<Real> &b)
    {
        if (m_full_rank) {
            mp_qrp->solve(x, b);
        } else {
#ifdef LS_WARN_RANK_DEFICIENT
            cout << "ls: A(" << m_rows << "x" << m_cols << ") is "
                 << "rank-deficient (rank = " << mp_svd->rank() 
                 << "), using SVD... "  << endl;
#endif
            mp_qrp->qTransposeB(b);
            mp_svd->solve(x, b);
            mp_qrp->pivot(x);
        }
    }
    
    void solve(double* const p_x, double* const p_b)
    {
        if (m_full_rank)
            mp_qrp->solve(p_x, p_b);
        else {
#ifdef LS_WARN_RANK_DEFICIENT
            cout << "ls: A(" << m_rows << "x" << m_cols << ") is "
                 << "rank-deficient (rank = " << mp_svd->rank() 
                 << "), using SVD... "  << endl;
#endif
            mp_qrp->qTransposeB(p_b);
            mp_svd->solve(p_x, p_b);
            mp_qrp->pivot(p_x);
        }
    }
    
    /**
     * Solves the square system A'Ax = b using the PR'RP'x = b where R and P
     * are the upper triangular and permutation matrices respectively resulting
     * from the QR decomposition of A with column pivoting.  Note that the rank
     * of A does not matter in this case.
     */
    void solveATA(Vector<Real> &x, Vector<Real> &b) {
        mp_qrp->solveATA(x, b);
    }
    
private:

    const int m_rows;
    const int m_cols;

    QRP<Real> *mp_qrp;
    SVD<Real> *mp_svd;
    
    bool m_full_rank;
};

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_LEASTSQUARES_H
