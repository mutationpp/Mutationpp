#ifndef PROBLEM_H
#define PROBLEM_H

#include <boost/numeric/odeint.hpp>
#include "mutation++.h"

#include "Data.h"

typedef boost::numeric::ublas::vector<double> vector_type;
using namespace boost::numeric::odeint;

class Problem {
public:
    Problem(Mutation::Mixture& l_mix, Data& l_data)
            : m_mix(l_mix),
              m_data(l_data){}
    virtual ~Problem(){}

    virtual vector_type computedxdt(const vector_type& x, const double t) = 0;
    int nEquations(){return m_data.nEquations();};

protected:
    Mutation::Mixture& m_mix;
    Data& m_data;
};

#endif /* PROBLEM_H */
