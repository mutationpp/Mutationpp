/**
 * @file NewtonSolver.h
 *
 * @brief Provides framework for a generalized Newton's method for solving
 * nonlinear equations.
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

#ifndef _NUMERICS_NEWTON_SOLVER_
#define _NUMERICS_NEWTON_SOLVER_

#include <cassert>
#include <typeinfo>
#include <iostream>

namespace Mutation {
    namespace Numerics {

//==============================================================================

/**
 * Implements the skeleton framework for Newton's method for solving a nonlinear
 * system of equations.  Classes should extend this class and provide the 
 * methods for computing f(x), J(x) = df/dx, inv(J)*f, and norm(f).
 */
template <typename T, typename Solver>
class NewtonSolver
{
public:

    /**
     * Constructor.
     */
    NewtonSolver();

    /**
     * Destructor.
     */
    virtual ~NewtonSolver() {}

    /**
     * Uses Newton's method to compute the zero of f(x) given an initial guess.
     */
    T& solve(T& x);

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
        assert(lag > 0);
        m_jacobian_lag = lag;
    }

    /**
     * Sets whether or not to write convergence history to the standard out.
     */
    void setWriteConvergenceHistory(bool hist) {
        m_conv_hist = hist;
    }

private:

    unsigned int m_max_iter;
    unsigned int m_jacobian_lag;
    double       m_epsilon;
    bool         m_conv_hist;

};

//==============================================================================

template <typename T, typename Solver>
NewtonSolver<T, Solver>::NewtonSolver()
    : m_max_iter(20),
      m_jacobian_lag(1),
      m_epsilon(1.0e-8),
      m_conv_hist(false)
{ }

//==============================================================================

template <typename T, typename Solver>
T& NewtonSolver<T, Solver>::solve(T& x)
{
    using std::cout;
    using std::endl;

    // Make sure jacobian is updated on first iteration
    unsigned int jac = m_jacobian_lag;

    // Compute the initial function value and its norm
    static_cast<Solver&>(*this).updateFunction(x);
    double f0_norm = 1.;// f0_norm = static_cast<Solver&>(*this).norm(); // 1.
    double resnorm = 1.0;

    // Let user know Newton's method is being called and print initial residual
    if (m_conv_hist) {
        cout << "Newton (" << typeid(Solver).name() << "): norm(f0) = "
             << f0_norm << endl;
    }

    // Iterate until converged
    for (int i = 0; resnorm > m_epsilon && i < m_max_iter; ++i) {
        if (m_conv_hist)
            cout << "  iter = " << i + 1;

        // Update jacobian if needed
        if (jac++ == m_jacobian_lag) {
            if (m_conv_hist)
                cout << ", update J";
            static_cast<Solver&>(*this).updateJacobian(x);
            jac = 1;
        }

        // x -= inv(J)*f
        x -= static_cast<Solver&>(*this).systemSolution();

        // Recompute f and norm(f)
        static_cast<Solver&>(*this).updateFunction(x);
        resnorm = static_cast<Solver&>(*this).norm() / f0_norm;

        if (m_conv_hist)
            cout << ", relative residual = " << resnorm << endl;
    }

    if (resnorm > m_epsilon && m_conv_hist) {
        cout << "Newton failed to converge after " << m_max_iter
             << " iterations with a relative residual of " << resnorm << endl;
    }
    return x;
}

//==============================================================================

    } // namespace Numerics
} // namespace Mutation

#endif // _NUMERICS_NEWTON_SOLVER_
