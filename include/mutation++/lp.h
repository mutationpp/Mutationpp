/**
 * @file lp.h
 *
 * @brief Implements the lp function for solving linear programming problems.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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

/**
 * In addition the enums LpResult and LpObjective are defined along with helper
 * routines simplex and simp1, simp2, and simp3 which are not meant to be used
 * as stand alone functions.
 *
 * @internal A test LP problem is solved by hand at the bottom of this file
 * which was used during testing of the simplex method.  Note that the simplex
 * function actually uses a more economic storage methodology for maintaining
 * the tableau than is represented in the test problem but it serves as a useful
 * case for understanding the basic algorithm.
 */
 
#ifndef NUMERICS_LP_H
#define NUMERICS_LP_H

#include <cassert>
#include <cmath>


#include <iostream>
#include <iomanip>

namespace Mutation {
    namespace Numerics {

/**
 * Enumerates possible outcomes of the lp function for solving linear
 * programming (LP) problems .
 *
 * @see lp()
 */
enum LpResult {
    SOLUTION_FOUND,     ///< lp() returns with a correct solution
    UNBOUNDED_SOLUTION, ///< the solution to the LP problem is unbounded
    NO_SOLUTION,        ///< there is no feasible solution to the LP problem
    LINEARLY_DEPENDENT  ///< solution exists but some contraints are lin. dep.
};

/**
 * Enumerates possible operations that the objective function can perform in a 
 * specified linear programming problem.
 *
 * @see lp()
 */
enum LpObjective {
    MINIMIZE,  ///< min_x( f'*x )
    MAXIMIZE   ///< max_x( f'*x )
};

template <typename Real>
int simplex(Real *const tableau, const int m, const int n, const int m1, 
            const int m2, int *const izrov, int *const iposv, const Real eps);


/**
 * Solves a linear programming problem.
 *
 * Solves a linear programming problem which finds x that either minimizes or
 * maximizes z = f'*x (where f is an N-dimensional vector and f' denotes the 
 * tranpose of f) subject to the constraints
 *
 *   A1*x <= b1 for A1 a M1 x N matrix
 *   A2*x >= b2 for A2 a M2 x N matrix
 *   A3*x  = b3 for A3 a M3 x N matrix
 *
 * and
 *
 *   x >= 0.
 *
 * The solution is obtained using the Simplex Algorithm.
 *
 * @param f          N dimensional objective function vector
 * @param objective  specifies which objective function to use (either max or 
 *                   min)
 * @param A          (M1+M2+M3) x N dimensional constraint matrix A = 
 *                   [A1' A2' A3']'
 * @param b          (M1+M2+M3) dimensional constraint vector b = [b1' b2' b3']'
 * @param m1         number of <= constraints
 * @param m2         number of >= constraints
 * @param x          N dimenional solution vector
 * @param z          solution of objective function
 *
 * @return Returns the possible outcomes of a linear programming problem
 *
 * @see LpObjective
 * @see LpResult
 * @see simplex()
 *
 * <b>Example usage:</b>
 * @code
 * // Solve the following LP problem:
 * // maximize z = x1 + x2 + 3*x3 - 0.5*x4 subject to
 * // x1 + 2*x3 <= 740
 * // 2*x2 - 7*x4 <= 0
 * // x2 - x3 + 2*x4 >= 0.5
 * // x1 + x2 + x3 + x4 = 9
 * // and x >= 0
 * // First set the problem up...
 * Matrix<double> A(4,4,0.0);
 * Vector<double> f(4,1.0);
 * Vector<double> b(4,0.0);
 * Vector<double> x(4);
 *  
 * A(0,0) = 1; A(0,2) = 2;
 * A(1,1) = 2; A(1,3) = -7;
 * A(2,1) = 1; A(2,2) = -1; A(2,3) = 2;
 * A(3,0) = 1; A(3,1) = 1; A(3,2) = 1; A(3,3) = 1;
 *  
 * f(2) = 3;
 * f(3) = -0.5;
 *  
 * b(0) = 740; 
 * b(2) = 0.5;
 * b(3) = 9;
 * 
 * // Now solve the problem...
 * // Upon return, x = [0.00, 3.33, 4.73, 0.95]' and z = 17.02
 * double z;
 * LpResult result = lp(f, MAXIMIZE, A, b, 2, 1, x, z);
 * @endcode
 *
 * @todo This function is probably better suited as a LpProblem class where the
 * initial tableau setup is completed in the constructor and a solve() method
 * could be used which actually performs the simplex method.  This way, it might
 * be possible to save large parts of the initial tableau for problems where 
 * the constraint vectors change.  Later if a sparse lp solver is added, the 
 * constructor for an LpProblem class could make the determination of which
 * solver to use based on the sparsity of the constraint matrix.
 *
 * @author J.B. Scoggins (jbscoggi@gmail.com)
 * @date   November 27, 2011
 */ 
//template <typename Real>
//LpResult lp(const Vector<Real> &f, const LpObjective &objective,
//            const Matrix<Real> &A, const Vector<Real> &b, const int m1,
//            const int m2, Vector<Real> &x, Real &z, int** iposv_ret = NULL)
//{
//    assert( f.size() == x.size() );
//    assert( A.rows() == b.size() );
//    assert( A.cols() == x.size() );
//
//    const int m = A.rows();
//    const int n = A.cols();
//
//    // Build the tableau for input into the simplex algorithm
//    Real tableau[m+2][n+1];
//    tableau[0][0] = static_cast<Real>(0);
//
//    // Fill tableau with f vector
//    for (int j = 1; j < n+1; ++j)
//        tableau[0][j] = (objective == MINIMIZE ? -f(j-1) : f(j-1));
//
//    // Fill tableau with A and b while ensuring that negative b values are
//    // handled according to the constraint type (note that the simplex method
//    // only works for positive b values)
//    int offset = 0;
//    Real temp;
//    for (int i = 1; i < m+1; ++i) {
//        tableau[i][0] = b(i-1);
//        // If b value for this row is negative
//        if (tableau[i][0] < static_cast<Real>(0)) {
//            // Change <= constraint to >= constraint
//            if (i <= m1) {
//                tableau[m1][0] = -tableau[i][0];
//                for (int j = 1; j < n+1; ++j)
//                    tableau[m1][j] = A(i-1,j-1);
//                offset++;
//            // Change >= constraint to <= constraint
//            } else if (i <= m1+m2) {
//                offset--;
//                // Prefer swap rather than have to call b(i-1) again
//                temp = -tableau[i][0];
//                tableau[i][0] = tableau[m1-offset][i];
//                tableau[m1-offset][0] = temp;
//                // Negate row of A
//                for (int j = 1; j < n+1; ++j) {
//                    tableau[i][j] = tableau[m1-offset][j];
//                    tableau[m1-offset][j] = A(i-1,j-1);
//                }
//            // Change sign of = constraint row
//            } else {
//                tableau[i][0] = -tableau[i][0];
//                for (int j = 1; j < n+1; ++j)
//                    tableau[i][j] = A(i-1,j-1);
//            }
//        // b value is positive so just copy as normal
//        } else {
//            // If this is a <= constraint then we have to keep in mind the
//            // offset before copying
//            if (i <= m1) {
//                for (int j = 1; j < n+1; ++j)
//                    tableau[i-offset][j] = -A(i-1,j-1);
//            // Otherwise we can copy as normal
//            } else {
//                for (int j = 1; j < n+1; ++j)
//                    tableau[i][j] = -A(i-1,j-1);
//            }
//        }
//    }
//
//    // Set the bottom row of the tableau to zero
//    for (int j = 0; j < n+1; ++j)
//        tableau[m+1][j] = static_cast<Real>(0);
//
//    // Compute a scaled epsilon based on the average b value
//    Real eps = static_cast<Real>(0);
//    for (int i = 1; i < m+1; ++i)
//        eps += tableau[i][0];
//
//    if (eps == static_cast<Real>(0))
//        eps = static_cast<Real>(1.0e-9);
//    else
//        eps *= static_cast<Real>(1.0e-9) / static_cast<Real>(m);
//
//    // Perform the simplex algorithm on the tableau
//    int izrov[n];
//    int* iposv = new int [m];
//    int icase = simplex(&tableau[0][0], m, n, m1-offset, m2+offset,
//                        izrov, iposv, eps);
//
//    //cout << "final tableau" << endl;
//    //for (int i = 0; i < m+2; ++i) {
//    //    for (int j = 0; j < n+1; ++j)
//    //        cout << setw(10) << tableau[i][j];
//    //    cout << endl;
//    //}
//    //cout << endl;
//
//    if (icase != 0) {
//        if (icase > 0)
//            return UNBOUNDED_SOLUTION;
//        else
//            return NO_SOLUTION;
//    }
//
//    // Now unravel the solution from the modified tableau
//    x = static_cast<Real>(0);
//    bool linind = false;
//    for (int i = 0; i < m; ++i) {
//        if (iposv[i] < n)
//            x(iposv[i]) = tableau[i+1][0];
//        else
//            linind = true;
//    }
//
//    if (iposv_ret != NULL)
//        *iposv_ret = iposv;
//    else
//        delete [] iposv;
//
//    // Also return value of z = min/max(f'*x)
//    z = tableau[0][0];
//
//    return (linind ? LINEARLY_DEPENDENT : SOLUTION_FOUND);
//} // lp


template <typename Real>
void simp1(const Real *const tableau, const int n, const int mm, 
           const int *const ll, const int nll, const bool use_abs, 
           int &kp, Real &bmax);
           
template <typename Real>
void simp2(const Real *const tableau, const int m, const int n, int &ip, 
           const int kp, const Real eps);
           
template <typename Real>
void simp3(Real *const tableau, const int n, const int i1, const int k1, 
           const int ip, const int kp);


/**
 * Performs the Simplex Algorithm on a prearranged tableau representing a linear
 * programming problem.  For a problem
 *
 *   max (f' * x)
 *    x
 *
 * subject to the constraints
 *
 *   Ax <= a for A M1 x N
 *   Bx >= b for B M2 x N
 *   Cx =  c for C M3 x N
 *
 * and
 *
 *    x >= 0, a >= 0, b >= 0, and c >= 0
 *
 * the tableau should be constructed as an M+2 x N+1 multidimensional array with
 * the entries
 *
 *  0   f'
 *  a  -A  
 *  b  -B
 *  c  -C
 *  0  0
 *
 * The code below is based on the Numerical Recipes (NR) routine simplx with 
 * minor modifications to remove the evil goto command and to adjust for the 
 * index shift between fortran and c/c++.  See the NR Chapter 10.8 for a
 * complete description of the simplx method.
 *
 * @param tableau  M+2xN+1-dim array filled with tableau data, on return holds
 *                 solution values in the first column in the order given by
 *                 iposv
 * @param m        total number of constraints
 * @param n        size of f and x vectors
 * @param m1       number of <= type constraints
 * @param m2       number of >= type constraints (m3 = m - m1 - m2)
 * @param izrov    N-dim array storing the remaining right-hand variable order
 *                 on return
 * @param iposv    M-dim array storing the left-hand variable order on return,
 *                 i.e. x[iposv[i]] = tableau[i+1][0] if iposv[i] < n+1
 *
 * @return  0 = success, < 0 = no solution, > 0 = solution unbounded
 *
 * @author J.B. Scoggins (jbscoggi@gmail.com)
 * @date   November 27, 2011
 */
template <typename Real>
int simplex(Real *const tableau, const int m, const int n, const int m1, 
            const int m2, int *const izrov, int *const iposv, const Real eps)
{
    const int m3   = m - m1 - m2;
    
    // Declare other variables
    Real bmax;
    int l1[n], nl1, kp, ip, i, j, is;
    bool need_exchange;
    
    // Initialize the list of columns admissible for exchange and make all
    // variables right-hand
    nl1 = n;
    for (i = 0; i < n; ++i) {
        l1[i]    = i;
        izrov[i] = i;
    }  
    
    // Initial left hand variables
    for (i = 0; i < m; ++i)
        iposv[i] = n+i;
    
    // PHASE 1:
    // Use auxiliary objective function to compute initial solution, if there
    // are no >= or = constraints then the origin is an initial solution
    if (m2 + m3 != 0) {
        // Initialize the list of m2 constraints whose slack variables have
        // never been exchanged out of the initial basis
        int l3[m2], kh;
        for (i = 0; i < m2; ++i)
            l3[i] = 1;
        
        // Compute the auxiliary objective function
        Real* a = &tableau[(m+1)*(n+1)];  // a points to row (m+1)
        std::fill(a, a+(n+1), static_cast<Real>(0));
        for (i = m1+1; i < m+1; ++i)
            for (j = 0; j < n+1; ++j)
                a[j] -= tableau[i*(n+1)+j];
        
        // Use simplex algorithm on auxiliary objective function
        while (true) {
            simp1(tableau, n, m+1, l1, nl1, false, kp, bmax);       

            // If auxiliary function is negative and can't be improved then no
            // feasible solution exists
            if (bmax <= eps && a[0] <= -eps)
                return -1;
                
            // If auxiliary function is zero and can't be improved then we have
            // a feasible starting vector
            if (bmax <= eps && a[0] <= eps) {
                need_exchange = false;
                for (ip = m1+m2; ip < m; ++ip) {
                    if (iposv[ip] == n+ip) {
                        simp1(tableau, n, ip, l1, nl1, true, kp, bmax);
                        if (bmax > eps) {
                            need_exchange = true;
                            break;
                        }
                    }
                }                

                if (!need_exchange) {
                    // Change sign of row for any m2 constraints still present
                    // from the initial basis
                    a = &tableau[(n+1)*(m1+1)];
                    for (i = m1; i < m1+m2; ++i, a += (n+1)) {
                        if (l3[i-m1] == 1) {
                            for (j = 0; j < n+1; ++j)
                                a[j] = -a[j];
                        }
                    }
                    
                    // Go to PHASE 2
                    break;
                }
            } else {
                // Locate the pivot element row
                simp2(tableau, m, n, ip, kp, eps);
                
                // Check to see if the auxiliary objective function is 
                // unbounded, if so then no feasible solution exists
                if (ip == -1)
                    return -1;
            }
            
            // Exchange a left- and a right-hand variable (phase one), then 
            // update the lists
            simp3(tableau, n, m+1, n, ip, kp);
            
            if (iposv[ip] >= n + m1 + m2) {
                for (j = 0; j < nl1; ++j) {
                    if (l1[j] == kp) break;
                }
                
                nl1--;
                for (is = j; is < nl1; ++is)
                    l1[is] = l1[is+1];
            } else {
                kh = iposv[ip] - m1 - n;
                if (kh >= 0 && l3[kh] != 0) {
                    l3[kh] = 0;
                    a[kp+1]++;
                    a = &tableau[0];
                    for (i = 0; i < m+2; ++i, a += (n+1))
                        a[kp+1] = -a[kp+1];
                    a -= (n+1);  // reset a to point to row (m+1)
                }
            }
            
            is = izrov[kp];
            izrov[kp] = iposv[ip];
            iposv[ip] = is;
        } // while (true)
    } // if (m2 + m3 != 0)
    
    // PHASE 2:
    // Found initial feasible solution, now optimize it
    while (true) {
        // Test z-row for doneness
        simp1(tableau, n, 0, l1, nl1, false, kp, bmax);
        
        if (bmax <= eps)
            return 0;
        
        // Locate a pivot element (phase two)
        simp2(tableau, m, n, ip, kp, eps);
        
        // Check for unbounded objective function
        if (ip == -1)
            return 1;
        
        // Exchange a left- and a right-hand variable (phase two) and update
        // lists
        simp3(tableau, n, m, n, ip, kp);
        
        is = izrov[kp];
        izrov[kp] = iposv[ip];
        iposv[ip] = is;
    }
    
} // simplex

/**
 * simplex auxiliary function 1.  Determines the maximum of those elements whose
 * index is contained in the supplied list ll, either with or without taking the
 * absolute value, as flagged by use_abs.
 */
template <typename Real>
void simp1(const Real *const tableau, const int n, const int mm, 
           const int *const ll, const int nll, const bool use_abs, 
           int &kp, Real &bmax)
{    
    // Only need to access row mm of the tableau
    const Real *const a = &tableau[mm*(n+1)];
    Real test;
    
    if (nll <= 0) {
        bmax = static_cast<Real>(0);
    } else {
        kp   = ll[0];
        bmax = a[kp+1];
        
        for (int k = 1; k < nll; ++k) {
            if (use_abs)
                test = std::abs(a[ll[k]+1]) - std::abs(bmax);
            else
                test = a[ll[k]+1] - bmax;
            
            if (test > 0) {            
                kp   = ll[k];
                bmax = a[kp+1];
            }                
        }        
    }
} // simp1

/**
 * simplex auxiliary function 2. Locate a pivot element, taking degeneracy into
 * account.
 */
template <typename Real>
void simp2(const Real *const tableau, const int m, const int n, int &ip, 
           const int kp, const Real eps)
{
    Real q, q0, q1, qp;
    
    int i;
    const Real* a = &tableau[(n+1)];
    for (i = 0; i < m; ++i, a += (n+1))
        if (a[kp+1] < -eps) break;
    
    // If there are no possible pivots, return
    if (i == m) {
        ip = -1;
        return;
    }
    
    q1 = -a[0] / a[kp+1];
    ip = i;
    
    a = &tableau[(n+1)*(ip+2)];
    for (i = ip+1; i < m; ++i, a += (n+1)) {
        if (a[kp+1] < -eps) {
            q = -a[0]/a[kp+1];
            if (q < q1) {
                ip = i;
                q1 = q;
            } else if (q == q1) {
                const Real* ap = &tableau[(n+1)*(ip+1)];
                for (int j = 0; j < n; ++j) {
                    qp = -ap[j+1] / ap[kp+1];
                    q0 = -a[j+1] / a[kp+1];
                    if (q0 != qp) break;
                }
                if (q0 < qp) ip = i;
            }
        }
    }
} // simp2

/**
 * simplex auxiliary function 3. Exchange a left and right hand variable.
 */
template <typename Real>
void simp3(Real *const tableau, const int n, const int i1, const int k1, 
           const int ip, const int kp)
{
    Real piv = static_cast<Real>(1) / tableau[(ip+1)*(n+1)+kp+1];

    // Update all elements except pivot row
    Real* a = &tableau[0];
    for (int i = 0; i < i1+1; ++i, a += (n+1)) {
        if (i != ip+1) {
            // Pivot column
            a[kp+1] *= piv;
            
            // Non-pivot column
            for (int j = 0; j < k1+1; ++j) {
                if (j != kp+1)
                    a[j] -= tableau[(ip+1)*(n+1)+j] * a[kp+1];
            }
        }
    }
    
    // Update the pivot row
    a = &tableau[(ip+1)*(n+1)];
    for (int j = 0; j < k1+1; ++j) {
        if (j != kp+1)
            a[j] *= -piv;
    }
    
    // Update the pivot element
    a[kp+1] = piv;
} // simp3
 
    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_LP_H

/*

The following list of tableau updates corresponds to the simplex method used on
the following linear programming problem:

maximize z = x1 + x2 + 3*x3 - 0.5*x4 for xi >= 0 and

x1 +      2*x3        <= 740
    2*x2 -       7*x4 <= 0
      x2 -  x3 + 2*x4 >= 0.5
x1 +  x2 +  x3 +   x4  = 9


Initial tableau:
================

                    x1       x2       x3       x4       y1       y2       y3
    ------------------------------------------------------------------------
 z  |    0.00     1.00     1.00     3.00    -0.50     0.00     0.00     0.00
z1  |  740.00    -1.00     0.00    -2.00     0.00    -1.00     0.00     0.00
z2  |    0.00     0.00    -2.00     0.00     7.00     0.00    -1.00     0.00
z3  |    0.50     0.00    -1.00     1.00    -2.00     0.00     0.00     1.00
z4  |    9.00    -1.00    -1.00    -1.00    -1.00     0.00     0.00     0.00
z'  | -749.50     2.00     4.00     2.00    -4.00     1.00     1.00    -1.00


PHASE 1:
========

                    x1       z2       x3       x4       y1       y2       y3
    ------------------------------------------------------------------------
 z  |    0.00     1.00    -0.50     3.00     3.00     0.00    -0.50     0.00
z1  |  740.00    -1.00    -0.00    -2.00     0.00    -1.00     0.00     0.00
x2  |    0.00     0.00    -0.50     0.00     3.50     0.00    -0.50     0.00
z3  |    0.50     0.00     0.50     1.00    -5.50     0.00     0.50     1.00
z4  |    9.00    -1.00     0.50    -1.00    -4.50     0.00     0.50     0.00
z'  | -749.50     2.00    -2.00     2.00    10.00     1.00    -1.00    -1.00

                    x1       z2       x3       z3       y1       y2       y3
    ------------------------------------------------------------------------
 z  |    0.27     1.00    -0.23     3.55    -0.55     0.00    -0.23     0.55
z1  |  740.00    -1.00     0.00    -2.00    -0.00    -1.00     0.00     0.00
x2  |    0.32     0.00    -0.18     0.64    -0.64     0.00    -0.18     0.64
x4  |    0.09     0.00     0.09     0.18    -0.18     0.00     0.09     0.18
z4  |    8.59    -1.00     0.09    -1.82     0.82     0.00     0.09    -0.82
z'  | -748.59     2.00    -1.09     3.82    -1.82     1.00    -0.09     0.82

                    x1       z2       z4       z3       y1       y2       y3
    ------------------------------------------------------------------------
 z  |   17.02    -0.95    -0.05    -1.95     1.05     0.00    -0.05    -1.05
z1  |  730.55     0.10    -0.10     1.10    -0.90    -1.00    -0.10     0.90
x2  |    3.32    -0.35    -0.15    -0.35    -0.35     0.00    -0.15     0.35
x4  |    0.95    -0.10     0.10    -0.10    -0.10     0.00     0.10     0.10
x3  |    4.72    -0.55     0.05    -0.55     0.45     0.00     0.05    -0.45
z'  | -730.55    -0.10    -0.90    -2.10    -0.10     1.00     0.10    -0.90

                    x1       z2       z4       z3       z1       y2       y3
    ------------------------------------------------------------------------
 z  |   17.02    -0.95    -0.05    -1.95     1.05    -0.00    -0.05    -1.05
y1  |  730.55     0.10    -0.10     1.10    -0.90    -1.00    -0.10     0.90
x2  |    3.32    -0.35    -0.15    -0.35    -0.35    -0.00    -0.15     0.35
x4  |    0.95    -0.10     0.10    -0.10    -0.10    -0.00     0.10     0.10
x3  |    4.72    -0.55     0.05    -0.55     0.45    -0.00     0.05    -0.45
z'  |   -0.00     0.00    -1.00    -1.00    -1.00    -1.00     0.00     0.00


PHASE 2:
========

                    x1       y2       y3
    ------------------------------------
 z  |   17.02    -0.95    -0.05    -1.05
y1  |  730.55     0.10    -0.10     0.90
x2  |    3.32    -0.35    -0.15     0.35
x4  |    0.95    -0.10     0.10     0.10
x3  |    4.72    -0.55     0.05    -0.45
z'  |   -0.00     0.00     0.00     0.00

z cannot be maximized any more so the solution is given as

x1 = 0.00
x2 = 3.32
x3 = 4.72
x4 = 0.95

and the max(z) = 17.02

*/
