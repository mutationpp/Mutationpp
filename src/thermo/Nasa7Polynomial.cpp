/**
 * @file Nasa7Polynomial.cpp
 *
 * @brief Implementation of Nasa7Polynomial class.
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

#include <cmath>

#include "Nasa7Polynomial.h"
#include "Utilities.h"

using namespace std;

namespace Mutation {
    namespace Thermodynamics {

Nasa7Polynomial::Nasa7Polynomial(const Nasa7Polynomial &to_copy)
{
    *this = to_copy;
}

Nasa7Polynomial &Nasa7Polynomial::operator=(const Nasa7Polynomial &to_copy)
{
    // Avoid self-assignment
    if (&to_copy != this) {    
        for (int i = 0; i < m_nr+1; i++)
            mp_tbounds[i] = to_copy.mp_tbounds[i];

        for (int i = 0; i < m_nr; ++i)
            for (int j = 0; j < m_nc; ++j)
                mp_coefficients[i][j] = to_copy.mp_coefficients[i][j];
    }
    
    // Enable cascading assignment
    return *this;
}

void Nasa7Polynomial::cp(const double *const p_params, double &cp) const
{
    int index = tRange(p_params[0]);
    
    cp = mp_coefficients[index][0] +
         mp_coefficients[index][1] * p_params[0] +
         mp_coefficients[index][2] * p_params[1] +
         mp_coefficients[index][3] * p_params[2] +
         mp_coefficients[index][4] * p_params[3];
}

void Nasa7Polynomial::enthalpy(const double *const p_params, double &h) const
{
    int index = tRange(p_params[0]);
    
    h = mp_coefficients[index][0] +
        mp_coefficients[index][1] * p_params[1] +
        mp_coefficients[index][2] * p_params[2] +
        mp_coefficients[index][3] * p_params[3] +
        mp_coefficients[index][4] * p_params[4] +
        mp_coefficients[index][5] * p_params[5];
}

void Nasa7Polynomial::entropy(const double *const p_params, double &s) const
{
    int index = tRange(p_params[1]);
    
    s = mp_coefficients[index][0] * p_params[0] +
        mp_coefficients[index][1] * p_params[1] +
        mp_coefficients[index][2] * p_params[2] +
        mp_coefficients[index][3] * p_params[3] +
        mp_coefficients[index][4] * p_params[4] +
        mp_coefficients[index][6];
}

void Nasa7Polynomial::gibbs(const double *const p_params, double &g) const
{
    int index = tRange(p_params[0]);
    
    g = mp_coefficients[index][0] * p_params[1] +
        mp_coefficients[index][1] * p_params[2] +
        mp_coefficients[index][2] * p_params[3] +
        mp_coefficients[index][3] * p_params[4] +
        mp_coefficients[index][4] * p_params[5] +
        mp_coefficients[index][5] * p_params[6] -
        mp_coefficients[index][6];
}

void Nasa7Polynomial::computeParams(const double &T, double *const params, 
    const ThermoFunction func)
{
    double T2, T3, T4;
    
    T2 = T*T;
    T3 = T2*T;
    T4 = T3*T;
    
    if (func == CP) {
        params[0] = T;
        params[1] = T2;
        params[2] = T3;
        params[3] = T4;
    }
    
    if (func == ENTHALPY) {
        params[0] = T;
        params[1] = T*0.5;
        params[2] = T2/3.0;
        params[3] = T3*0.25;
        params[4] = T4*0.2;
        params[5] = 1.0/T;
    }
    
    if (func == ENTROPY) {
        params[0] = log(T);
        params[1] = T;
        params[2] = T2*0.5;
        params[3] = T3/3.0;
        params[4] = T4*0.25;
    }
    
    if (func == GIBBS) {
        params[0] = T;
        params[1] = 1.0-log(T);
        params[2] = -T/2.0;
        params[3] = -T2/6.0;
        params[4] = -T3/12.0;
        params[5] = -T4/20.0;
        params[6] = 1.0/T;
    }
}

std::istream &operator>>(istream &in, Nasa7Polynomial &n7)
{
    // Read the first line and extract the temperature limits
    char temp[100];
    in.get(temp, 49);
    in.get(temp, 10); n7.mp_tbounds[0] = atof(temp);
    in.get(temp, 10); n7.mp_tbounds[2] = atof(temp);
    in.get(temp, 10); n7.mp_tbounds[1] = atof(temp);
    in.getline(temp, 100);
    
    // Second line
    for (int i = 0; i < 5; ++i) {
        in.get(temp, 16); 
        n7.mp_coefficients[1][i] = atof(temp);
    }
    in.getline(temp, 100);
    
    // Third line
    for (int i = 5; i < 7; ++i) {
        in.get(temp, 16); 
        n7.mp_coefficients[1][i] = atof(temp);
    }
    for (int i = 0; i < 3; ++i) {
        in.get(temp, 16); 
        n7.mp_coefficients[0][i] = atof(temp);
    }
    in.getline(temp, 100);
    
    // Fourth line
    for (int i = 3; i < 7; ++i) {
        in.get(temp, 16);
        n7.mp_coefficients[0][i] = atof(temp);
    }
    in.getline(temp, 100);
    
    return in;
}

    } // namespace Thermodynamics
} // namespace Mutation

