/**
 * @file Nasa9Polynomial.cpp
 *
 * @brief Implementation of Nasa9Polynomial class.
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
#include <cstdlib>
#include <string>
#include "Nasa9Polynomial.h"
#include "Utilities.h"

using namespace std;
using namespace Mutation::Utilities;

namespace Mutation {
namespace Thermodynamics {

Nasa9Polynomial::Nasa9Polynomial()
    : m_nr(0), m_nc(0), mp_coefficients(NULL), mp_tbounds(NULL)
{ }

Nasa9Polynomial::Nasa9Polynomial(const Nasa9Polynomial& to_copy)
    : m_nr(to_copy.m_nr), m_nc(to_copy.m_nc), 
      mp_coefficients(
        to_copy.mp_coefficients != NULL ? new double* [m_nr] : NULL), 
      mp_tbounds(
        to_copy.mp_tbounds != NULL ? new double [m_nr+1] : NULL)
{
    if (mp_coefficients != NULL) {
        for (int i = 0; i < m_nr; ++i) {
            mp_coefficients[i] = new double [m_nc];
            for (int j = 0; j < m_nc; ++j)
                mp_coefficients[i][j] = to_copy.mp_coefficients[i][j];
        }
    }
    
    if (mp_tbounds != NULL) {
        for (int i = 0; i < m_nr+1; ++i)
            mp_tbounds[i] = to_copy.mp_tbounds[i];
    }
}

Nasa9Polynomial::~Nasa9Polynomial()
{
    if (mp_coefficients != NULL) {
        for (int i = 0; i < m_nr; ++i)
            delete [] mp_coefficients[i];
        delete [] mp_coefficients;
    }
    
    if (mp_tbounds != NULL)
        delete [] mp_tbounds;
}

Nasa9Polynomial& Nasa9Polynomial::operator = (Nasa9Polynomial left)
{
    swap(*this, left);
    return *this;
}

void Nasa9Polynomial::cp(const double *const p_params, double &cp) const
{
    int tr = tRange(p_params[3]);

    cp = mp_coefficients[tr][0] * p_params[0];
    for (int i = 1; i < 7; ++i)
        cp += mp_coefficients[tr][i] * p_params[i];
}

void Nasa9Polynomial::enthalpy(const double *const p_params, double &h) const
{
    int tr = tRange(2.0 * p_params[3]);

    h = mp_coefficients[tr][0] * p_params[0];
    for (int i = 1; i < 8; ++i)
        h += mp_coefficients[tr][i] * p_params[i];
}

void Nasa9Polynomial::entropy(const double *const p_params, double &s) const
{
    int tr = tRange(p_params[3]);

    s = mp_coefficients[tr][8];
    for (int i = 0; i < 7; ++i)
        s += mp_coefficients[tr][i] * p_params[i];
}

void Nasa9Polynomial::gibbs(const double *const p_params, double &g) const
{
    int tr = tRange(-2.0 * p_params[3]);

    g = -mp_coefficients[tr][8];
    for (int i = 0; i < 8; ++i)
        g += mp_coefficients[tr][i] * p_params[i];
}

void Nasa9Polynomial::computeParams(const double &T, double *const params, 
    const ThermoFunction func)
{
    double T2, T3, T4;    
    
    T2 =  T * T;
    T3 = T2 * T;
    T4 = T3 * T;
    
    // Cp(T)/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
    if (func == CP) {
        params[0] = 1.0 / T2;
        params[1] = 1.0 / T;
        params[2] = 1.0;
        params[3] = T;
        params[4] = T2;
        params[5] = T3;
        params[6] = T4;
    }
    
    // H(T)/RT = -a1/T^2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + 
    //           a7*T^4/5 + b1/T
    if (func == ENTHALPY) {
        params[0] = -1.0 / T2;
        params[1] = log(T) / T;
        params[2] = 1.0;
        params[3] = 0.5 * T;
        params[4] = T2 / 3.0;
        params[5] = 0.25 * T3;
        params[6] = T4 / 5.0;
        params[7] = 1.0 / T;
    }
    
    // S(T)/R = –a1/(2*T^2) – a2/T + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3 + 
    //          a7*T^4/4 + b2
    if (func == ENTROPY) {
        params[0] = -0.5 / T2;
        params[1] = -1.0 / T;
        params[2] = log(T);
        params[3] = T;
        params[4] = 0.5 * T2;
        params[5] = T3 / 3.0;
        params[6] = 0.25 * T4;
    }
    
    // G(T)/RT = H(T)/RT - S(T)/R
    if (func == GIBBS) {
        params[0] = -0.5 / T2;
        params[1] = (log(T) + 1.0) / T;
        params[2] = 1.0 - log(T);
        params[3] = -0.5 * T;
        params[4] = -T2 / 6.0;
        params[5] = -T3 / 12.0;
        params[6] = -T4 / 20.0;
        params[7] = 1.0 / T;
    }
}

int Nasa9Polynomial::tRange(double T) const
{
    for (int i = 1; i < m_nr; ++i)
        if (T < mp_tbounds[i]) return (i-1);
    return (m_nr-1);
}

std::istream& operator >> (std::istream& in, Nasa9Polynomial& n9)
{
    string line;
    string name;
    
    // First line contains:
    // cols   format      value
    // 1-24   A24         Species name or formula
    // 25-80  A56         Comments and data source
    getline(in, line);
    name = line.substr(0, min(static_cast<int>(line.find(" ")), 24));
    
    // Second line contains:
    // cols   format      value
    // 1-2    I2          Number of T intervals
    // 4-9    A6          Optional identification code
    // 11-50  5(A2,F6.2)  Chemical formulasymbols and numbers (all capitals)
    // 51-52  I2          Zero for gas, nonzero for condensed
    // 53-65  F13.7       Molecular weight
    // 66-80  F15.5       Heat of formation at 298.15 K, J/mol
    getline(in, line);
    n9.m_nr = atoi(line.substr(0,2).c_str());
    
    // Allocate the storage (clean up storage first if already allocated)
    n9.m_nc = 9;
    
    if (n9.mp_tbounds != NULL)
        delete [] n9.mp_tbounds;
        
    n9.mp_tbounds = new double [n9.m_nr+1];
    
    if (n9.mp_coefficients != NULL) {
        for (int i = 0; i < n9.m_nr; ++i)
            delete [] n9.mp_coefficients[i];
        delete [] n9.mp_coefficients;
    }
    
    n9.mp_coefficients = new double* [n9.m_nr];
    for (int i = 0; i < n9.m_nr; ++i)
        n9.mp_coefficients[i] = new double [n9.m_nc];
    
    // Read in the polynomial for each temperature range
    for (int i = 0; i < n9.m_nr; ++i) {
        // Third line contains:
        // cols   format      value
        // 2-21   2F11.3      Temperature range
        // 23     I1          Number of coefficients for Cp(T)/R 
        // 24-63  8F5.1       T exponents in empirical equation for Cp(T)/R 
        // 66-80  F15.3       [H(298.15) – H(0)], J/mol
        getline(in, line);
        double t1 = atof(line.substr( 1,10).c_str()); // T_{1,i}
        double t2 = atof(line.substr(11,10).c_str()); // T_{2,i}
        
        if (t2 <= t1) {
            throw InvalidInputError("NASA-9 polynomial", name)
                << "In the " << i+1 << ordinalSuffix(i+1) << " temperature "
                << "range for species " << name << ", "
                << "T2 (" << t2 << ") is less than T1 (" << t1 << ").";
        }
        
        if (i == 0) {
            n9.mp_tbounds[i]   = t1;
            n9.mp_tbounds[i+1] = t2;
        } else {
            if (t1 != n9.mp_tbounds[i]) {
                throw InvalidInputError("NASA-9 polynomial", name)
                    << "T1 of " << i+1 << ordinalSuffix(i+1) << " range " << i << " does not match T2 of "
                    << "range " << i-1 << " for species " << name << ".";
            }
            
            n9.mp_tbounds[i+1] = t2;
        }
        
        // Verify that there are 7 exponents because that is an assumption we
        // are making
        if (atoi(line.substr(22,1).c_str()) != 7) {
            throw InvalidInputError("NASA-9 polynomial", name)
                << "Expecting 7 exponents in the polynomial.";
        }
        
        // Verify that the exponents are -2, -1, 0, 1, 2, 3, 4 because that is
        // an assumption that we make
        bool correct_exp = true;
        correct_exp &= (atof(line.substr(23,5).c_str()) == -2.0f);
        correct_exp &= (atof(line.substr(28,5).c_str()) == -1.0f);
        correct_exp &= (atof(line.substr(33,5).c_str()) ==  0.0f);
        correct_exp &= (atof(line.substr(38,5).c_str()) ==  1.0f);
        correct_exp &= (atof(line.substr(43,5).c_str()) ==  2.0f);
        correct_exp &= (atof(line.substr(48,5).c_str()) ==  3.0f);
        correct_exp &= (atof(line.substr(53,5).c_str()) ==  4.0f);
        
        if (!correct_exp) {
            throw InvalidInputError("NASA-9 polynomial", name)
                << "Only exponents (-2, -1, 0, 1, 2, 3, 4) are supported.";
        }
    
        // Load the coefficients from lines 4 and 5.  Note that there is a 
        // replace() in each line which converts the "D" in the string to an "E"
        // so that atof() will work properly.
        
        // Fourth line contains:
        // cols   format      value
        // 1-80   5D16.9      First five coefficients for Cp(T)/R
        getline(in, line);
        n9.mp_coefficients[i][0] = atof(line.substr( 0,16).replace(12,1,"E").c_str()); // a_{1,i}
        n9.mp_coefficients[i][1] = atof(line.substr(16,16).replace(12,1,"E").c_str()); // a_{2,i}
        n9.mp_coefficients[i][2] = atof(line.substr(32,16).replace(12,1,"E").c_str()); // a_{3,i}
        n9.mp_coefficients[i][3] = atof(line.substr(48,16).replace(12,1,"E").c_str()); // a_{4,i}
        n9.mp_coefficients[i][4] = atof(line.substr(64,16).replace(12,1,"E").c_str()); // a_{5,i}
        
        // Fifth line contains:
        // cols   format      value
        // 1-32   2D16.9      Last two coefficients for Cp(T)/R 
        // 49-80  2D16.9      Integration constants b1 and b2 
        getline(in, line);
        n9.mp_coefficients[i][5] = atof(line.substr( 0,16).replace(12,1,"E").c_str()); // a_{6,i}
        n9.mp_coefficients[i][6] = atof(line.substr(16,16).replace(12,1,"E").c_str()); // a_{7,i}
        n9.mp_coefficients[i][7] = atof(line.substr(48,16).replace(12,1,"E").c_str()); // b_{1,i}
        n9.mp_coefficients[i][8] = atof(line.substr(64,16).replace(12,1,"E").c_str()); // b_{2,i}
    }
    
//    // We have now loaded all the data we need, write out the data to check that
//    // it is working
//    cout << "Species: " << name << endl;
//    for (int i = 0; i < n9.m_nr; ++i) {
//        cout << "T-range: " << n9.mp_tbounds[i] << " - " << n9.mp_tbounds[i+1] << endl;
//        for (int j = 0; j < n9.m_nc; ++j)
//            cout << "a" << j << ": " << n9.mp_coefficients[i][j] << endl; 
//    }

    return in;
}

void swap(Nasa9Polynomial& left, Nasa9Polynomial& right)
{
    using std::swap;
    swap(left.m_nr, right.m_nr);
    swap(left.m_nc, right.m_nc);
    swap(left.mp_coefficients, right.mp_coefficients);
    swap(left.mp_tbounds, right.mp_tbounds);
}

} // namespace Thermodynamics
} // namespace Mutation

