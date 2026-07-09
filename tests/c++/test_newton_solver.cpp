/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include <iostream>

#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include "NewtonSolver.h"
#include <catch.hpp>
#include <Eigen/Dense>

//==============================================================================

/// Helper class for Newton solver tests
class NewtonSolverTest :
    public Mutation::Numerics::NewtonSolver<
        Eigen::VectorXd, NewtonSolverTest>
{
public:
    NewtonSolverTest()
        : mv_X(2),
          mv_dX(2),
          mv_f(2),
          mv_f_old(2),
          mv_df(2),
          mv_trans(1, 2),
          m_jac(2, 2),
          m_inv_jac(2, 2),
          m_tol(1.e-12)
    {
        // Setup NewtonSolver
        setMaxIterations(10);
        setWriteConvergenceHistory(true);
        setEpsilon(m_tol);
    }

    void updateFunction(Eigen::VectorXd& v_X) {
        // Save the old f
        mv_f_old = mv_f;

        // g(0) = 10.0*(x2 - x1^2)
        // g(1) = 1.0 - x1
        mv_f(0) = 10.0*(v_X(1) - v_X(0)*v_X(0));
        mv_f(1) = 1.0 - v_X(0);

        // Compute and store difference from new to old
        mv_df = mv_f - mv_f_old;
    }

    void updateJacobian(Eigen::VectorXd& v_X)
    {
        m_jac(0, 0) = -20.0*v_X(0);
        m_jac(0, 1) = 10.0;
        m_jac(1, 0) = -1.0;
        m_jac(1, 1) = 0.0;
    }

    Eigen::VectorXd& systemSolution()
    {
        mv_dX = m_jac.fullPivLu().solve(mv_f);
        return mv_dX;
    }

    double norm() {
        return mv_dX.lpNorm<Eigen::Infinity>();
    }

//==============================================================================

private:
    Eigen::VectorXd mv_X;
    Eigen::VectorXd mv_dX;
    Eigen::VectorXd mv_f;
    Eigen::VectorXd mv_f_old;
    Eigen::VectorXd mv_df;
    Eigen::MatrixXd mv_trans;
    Eigen::MatrixXd m_jac;
    Eigen::MatrixXd m_inv_jac;
    double m_tol;
};

//==============================================================================

/**
 * Checks that Newton's method solver works by solving the Rosenbrock function
 */
TEST_CASE("Newton solver returns the correct root",
    "[numerics][NewtonSolver]"
)
{
    constexpr double tol = 1.0e-12;
    NewtonSolverTest test;
    Eigen::VectorXd x(2);
    // Set the initial guess
    x(0) = 10.0;
    x(1) = 10.0;
    // Solve
    bool converged = test.solve(x);
    // Make sure it converged
    CHECK(converged);
    // Make sure the roots are 1.0
    INFO("x:\n" << x);
    CHECK(x(0) == Catch::Detail::Approx(1.0).margin(tol));
    CHECK(x(1) == Catch::Detail::Approx(1.0).margin(tol));
}

//==============================================================================

