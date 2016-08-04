#ifndef BASICOUTPUT_H
#define BASICOUTPUT_H

#include <iostream>
#include <boost/numeric/odeint.hpp>

#include "Data.h"

typedef boost::numeric::ublas::vector<double> vector_type;

class BasicOutput {
public:
    BasicOutput(const Data& l_data) // prints results at each step
           : n_eq(l_data.nEquations()),
             n_steps(0)
    {   steps_counter = 0;   }

    BasicOutput(const Data& l_data, size_t l_n_steps) // Prints results each n steps
           : n_eq(l_data.nEquations()),
             n_steps(l_n_steps)
    {   steps_counter = 0;   }

    ~BasicOutput(){}

    void operator()(const vector_type& l_sol, const double x){
       
        steps_counter++; 
        if(steps_counter >= n_steps) {

            steps_counter = 0;   // restart the counter
            std::cout << "At x = " << x << "  ";
            for (int i_eq = 0; i_eq < n_eq; i_eq++){
                std::cout << l_sol[i_eq] << "  ";
            }
            std::cout << std::endl;
        }
    }
private:
    const size_t n_eq;
    const size_t n_steps;
    size_t steps_counter;

};

#endif /* BASICOUTPUT_H */
