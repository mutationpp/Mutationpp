#ifndef JACOBIAN_H
#define JACOBIAN_H

typedef boost::numeric::ublas::matrix< double > matrix_type;

#include "Problem.h"

class jacobian {
public:
    jacobian(Problem& l_problem)
            : m_problem(l_problem),
              m_neq(l_problem.nEquations()),
              v_x(l_problem.nEquations()),
              dxdt(m_neq, 0.0),
              dxdt_pert(m_neq, 0.0),
              x_unpert(0.0),
              m_pert(1.e-6),
              m_tol(1.e-16){}

    void operator()(const vector_type& x, matrix_type& J, const double& t, vector_type& dfdt){

        dxdt = m_problem.computedxdt(x, t);
        v_x = x;

        for (int j = 0; j < m_neq; ++j){
            x_unpert = v_x(j);
            v_x(j) += v_x(j) * m_pert;
            dxdt_pert = m_problem.computedxdt(v_x, t);

            for (int i = 0; i < m_neq; ++i){
                J(i, j) = (dxdt_pert(i)-dxdt(i)) / ( x_unpert*m_pert + m_tol );
            }

            v_x(j) = x_unpert;
            dfdt[j] = 0.0;
        }

    }

    void setPertParam(const double& l_pert){m_pert = l_pert;}

private:
    Problem& m_problem;
    size_t m_neq;

    vector_type v_x;
    vector_type dxdt, dxdt_pert;
    double x_unpert;
    double m_pert, m_tol;
};

#endif /* JACOBIAN_H */
