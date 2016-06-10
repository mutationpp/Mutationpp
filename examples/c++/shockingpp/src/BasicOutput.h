#ifndef BASICOUTPUT_H
#define BASICOUTPUT_H

#include <iostream>
#include <boost/numeric/odeint.hpp>

#include "Data.h"

typedef boost::numeric::ublas::vector<double> vector_type;

class BasicOutput {
public:
    BasicOutput(const Data& l_data): n_eq(l_data.nEquations()){}
    ~BasicOutput(){}

    void operator()(const vector_type& l_sol, const double x){
        std::cout << "At x = " << x << "  ";
        for (int i_eq = 0; i_eq < n_eq; i_eq++){
            std::cout << l_sol[i_eq] << "  ";
        }
        std::cout << std::endl;
    }
private:
    const size_t n_eq;
};

#endif /* BASICOUTPUT_H */
