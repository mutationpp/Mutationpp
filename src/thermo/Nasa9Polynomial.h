/**
 * @file Nasa9Polynomial.h
 *
 * @brief Declaration of Nasa9Polynomial class.
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

#ifndef THERMO_NASA_9_POLYNOMIAL_H
#define THERMO_NASA_9_POLYNOMIAL_H

#include <iostream>

namespace Mutation {
    namespace Thermodynamics {

/**
 * Implements the NASA 9-coefficient polynomial model for a single species
 * thermodynamic properties.
 */
class Nasa9Polynomial
{
public:
    
    /**
     * Constructor, empty polynomial.
     */
    Nasa9Polynomial();
    
    /**
     * Copy constructor.
     */
    Nasa9Polynomial(const Nasa9Polynomial& to_copy);
    
    /**
     * Destructor, frees the polynomial data.
     */
    ~Nasa9Polynomial();
    
    /**
     * Assignment operator.
     */
    Nasa9Polynomial& operator = (Nasa9Polynomial left);
    
    static int nCoefficients() { return 9; }
    
    double minT() const {
        return mp_tbounds[0];
    }
    
    double maxT() const {
        return mp_tbounds[m_nr];
    }
    
    /**
     * Computes dimensionless specific heat Cp/Ru.
     * @see computeParams()
     */
    void cp(const double *const p_params, double &cp) const;
    
    /**
     * Computes dimensionless enthalpy H/Ru/T.
     * @see computeParams()
     */
    void enthalpy(const double *const p_params, double &h) const;
    
    /**
     * Computes dimensionless entropy S/Ru.
     * @see computeParams()
     */
    void entropy(const double *const p_params, double &s) const;
    
    /**
     * Computes dimensionless Gibbs free energy G/Ru/T.
     * @see computeParams()
     */
    void gibbs(const double *const p_params, double &g) const;
    
    /**
     * Enumerates functions that require a parameters list.
     * @see computeParams()
     */
    enum ThermoFunction
    {
        CP,
        ENTHALPY,
        ENTROPY,
        GIBBS,
        THERMO
    };
    
    /**
     * Computes the constant parameters as functions of temperature that must
     * be passed to the individual thermodynamic functions. The params array is
     * filled with the parameters needed by the function specified by func.  The
     * dimension of params should be at least as long as the following depending
     * on the value of func:
     *
     *  CP       - params[7]    \n
     *  ENTHALPY - params[8]    \n
     *  ENTROPY  - params[7]    \n
     *  GIBBS    - params[8]    \n
     *  THERMO   - params[22]
     *
     * Computing the parameters this way saves computation when they are reused
     * to compute thermodynamic properties for many species at once.
     *
     * @param T       the temperature
     * @param params  on return, equals the params list needed by the 
     *                appropriate thermodynamic function
     * @param func    specifices which thermodynamic function to use
     *
     * @see ThermoFunction
     */
    static void computeParams(const double &T, double *const params, 
                              const ThermoFunction func);
    
    friend std::istream& operator >> (std::istream& in, Nasa9Polynomial& n9);
    friend void swap(Nasa9Polynomial& left, Nasa9Polynomial& right);
    
private:

    /**
     * Determines which temperature range the temperature falls into.
     */
    int tRange(double T) const;

    int m_nr;
    int m_nc;

    double **mp_coefficients;
    double *mp_tbounds;
    
};

/**
 * Instantiates a Nasa9Polynomial object from an input stream with the position
 * starting exactly at the beginning of the first line of the polynomial data.
 */
std::istream& operator >> (std::istream& in, Nasa9Polynomial& n9); 

/**
 * Swaps one Nasa9Polynomial object with another.
 */
void swap(Nasa9Polynomial& left, Nasa9Polynomial& right);

} // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_NASA_9_POLYNOMIAL_H
