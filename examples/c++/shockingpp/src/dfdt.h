#ifndef DFDT_H
#define DFDT_H

class dfdt {
public:
    dfdt(Problem& l_problem)
            : m_problem(l_problem){} //@todo Initialize Problem

    void operator()(const vector_type& x, vector_type& dxdt, const double t){
        dxdt = m_problem.computedxdt(x, t);
    }

private:
    Problem& m_problem;
};

#endif /* DFDT_H */
