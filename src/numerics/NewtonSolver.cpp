//
//  NewtonSolver.cpp
//  Mutation++
//
//  Created by JB Scoggins on 8/29/13.
//  Copyright (c) 2013 JB Scoggins. All rights reserved.
//

#include "NewtonSolver.h"
#include <iostream>

using std::cout;
using std::endl;

namespace Mutation {
    namespace Numerics {
    
//==============================================================================

NewtonSolver::NewtonSolver(unsigned int n)
    : m_n(n),
      m_max_iter(20),
      m_jacobian_lag(1),
      m_epsilon(1.0e-8),
      m_conv_hist(false)
{ }

//==============================================================================

void NewtonSolver::solve(Vector<double>& x)
{
    Vector<double> f(m_n);
    Vector<double> dx(m_n);
    Matrix<double> J(m_n,m_n);
    
    unsigned int jac = m_jacobian_lag;
    
    // Compute the initial function value and its norm
    computeFunction(x, f);
    double f0_norm = f.norm2();
    double resnorm = 1.0;
    
    if (m_conv_hist) {
        cout << "Newton: norm(f0) = " << f0_norm << endl;
    }
    
    // Iterate until converged
    for (int i = 0; resnorm > m_epsilon && i < m_max_iter; ++i) {
        if (m_conv_hist)
            cout << "  iter = " << i + 1;
        
        // Update jacobian if needed
        if (jac++ == m_jacobian_lag) {
            if (m_conv_hist)
                cout << ", update J";
            computeJacobian(x, J);
            jac = 1;
        }
        
        // Solve the linear system J*dx = f
        systemSolution(dx, J, f);
        
        // Update the solution
        x -= dx;
        
        // Recompute f and norm(f)
        computeFunction(x, f);
        resnorm = f.norm2() / f0_norm;
        
        if (m_conv_hist)
            cout << ", relative residual = " << resnorm << endl;
    }
    
    if (resnorm > m_epsilon) {
        cout << "Newton failed to converge after " << m_max_iter
             << " iterations with a relative residual of " << resnorm << endl;
    }
}

//==============================================================================

    } // namespace Numerics
} // namespace Mutation
