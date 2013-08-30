//
//  NewtonSolver.h
//  Mutation++
//
//  Created by JB Scoggins on 8/29/13.
//  Copyright (c) 2013 JB Scoggins. All rights reserved.
//

#ifndef __Mutation____NewtonSolver__
#define __Mutation____NewtonSolver__

#include <cassert>
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {

/**
 * Implements the skeleton framework for Newton's method for solving a nonlinear
 * system of equations.  Classes should extend this class and provide the 
 * methods for computing f(x), J(x) = df/dx, and inv(J)*f.
 */
class NewtonSolver
{
public:
    
    /**
     * Constructor.
     */
    NewtonSolver(unsigned int n);
    
    /**
     * Destructor.
     */
    virtual ~NewtonSolver() {}
    
    /**
     * Uses Newton's method to compute the zero of f(x).
     */
    void solve(Vector<double>& x);
    
    /**
     * Set the residual norm tolerance.
     */
    void setEpsilon(const double eps) {
        assert(eps > 0.0);
        assert(eps < 1.0);
        m_epsilon = eps;
    }
    
    /**
     * Set the maximum number of iterations to use.
     */
    void setMaxIterations(const unsigned int iters) {
        m_max_iter = iters;
    }
    
    /**
     * Set the number of iterations to lag the Jacobian update.
     */
    void setJacobianLag(const unsigned int lag) {
        m_jacobian_lag = lag;
    }
    
    /**
     * Sets whether or not to write convergence history to the standard out.
     */
    void setWriteConvergenceHistory(bool hist) {
        m_conv_hist = hist;
    }

protected:

    /**
     * Should be implemented by the user class to provide the function
     * evaluation f(x).
     */
    virtual void computeFunction(
        const Vector<double>& x, Vector<double>& f) = 0;
    
    /**
     * Should be implemented by the user class to provide the Jacobian 
     * evaluation J = df/dx.
     */
    virtual void computeJacobian(
        const Vector<double>& x, Matrix<double>& f) = 0;
    
    /**
     * Should be implemented by the user class to solve the linear system
     * J*dx = f.
     */
    virtual void systemSolution(
        Vector<double>& dx, const Matrix<double>& J, Vector<double>& f) = 0;

private:

    unsigned int m_n;
    unsigned int m_max_iter;
    unsigned int m_jacobian_lag;
    
    double m_epsilon;
    
    bool m_conv_hist;
    
};

    } // namespace Numerics
} // namespace Mutation

#endif /* defined(__Mutation____NewtonSolver__) */
