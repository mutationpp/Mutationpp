/**
 * @file Numerics.h
 *
 * @brief Defines several useful numerical utilities, such as integration
 * 	  and factorials.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef NUMERICS
#define NUMERICS

#include <cmath>
#include <functional>

namespace Mutation {
    namespace Numerics {

    //! \brief integrate a function f on the semi-infinite interval [a,+infinity)
    //
    // it performs a change of variables and then integrates over the interval [0,1]
    // the interval is split into equal-sized subintervals and the integral is computed over each
    // subinterval, the amount of subintervals is determined by the subdivisions parameter
    // if error_estimate is not a NULL pointer, write the error_estimate to the passed variable
    //
    // default values: int subdivisions=5, double a = 0.0, double* error_estimate = nullptr
    //
    // usage: to integrate a function f(double x1, double x2, int x3)
    // over x1 from 0 to infinity (x2, x3 are fixed values defined somewhere in the code):
    //
    // auto integrand = [x2, x3](double x1) {return f(x1, x2, x3)};
    // double result = kappa::integrate_semi_inf(integrand);
    double integrate_semi_inf(std::function<double(double)>f, double a=0, int subdivisions=5, double* error_estimate=nullptr);

    // integrate a function f on the finite interval [a, b], where b > a
    // if error_estimate is not a NULL pointer, write the error_estimate to the passed variable
    //
    // default values: double* error_estimate = nullptr
    //
    // usage: to integrate a function f(double x1, double x2, int x3)
    // over x1 from -1 to 1 (x2, x3 are fixed values defined somewhere in the code):
    //
    // auto integrand = [x2, x3](double x1) {return f(x1, x2, x3)};
    // double result = kappa::integrate_semi_interval(integrand, -1, 1);
    double integrate_interval(std::function<double (double)>f, double a, double b, double* error_estimate=nullptr);

    // find the value of a parameter i of type int such that f(i) < max_value, f(i+1) > max_value
    // the search starts from i=start and then i is increased by 1 on each iteration (up to INT_MAX-2)
    // returns -1 if no value is found up to INT_MAX-2
    //
    // default values: int start=0
    int find_max_value(std::function<double(int)>f, double max_value, int start=0);

    // Calculate n!
    double factorial(int n);

    // Calculate end! / start!, where start and end are integers
    double fact_div_fact(int start, int end);

    // Convert a value given in cm^-1 to Joules
    double convert_cm_to_Joule(double x);

    // precomputed values of factorials from 0! to 69!, factorial_table[i] = i!
    //const extern arma::vec::fixed<70> factorial_table; 

    //arma::vec::fixed<70> compute_factorial_table(void);

    // coefficients for Born-Mayer potential, Kustova-Nagnibeda (5.100)
    constexpr double Born_Mayer_coeff_array[6][4] = { 	-267.0, 201.570, 174.672, 54.305,
                          				26700, -19226.5, -27693.8, -10860.9,
				                      	-8.9e5, 6.3201e5, 1.0227e6, 5.4304e5,
                          				-33.0838, 20.0862, 72.1059, 68.5001,
                          				101.571, -56.4472, -286.393, -315.4531,
                          				-87.7036, 46.3130, 277.146, 363.1807 }; 
    }
}

#endif 
