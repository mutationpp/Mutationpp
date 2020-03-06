/**
 * @file Nasa7Polynomial.h
 *
 * @brief Declaration of Nasa7Polynomial class.
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


#ifndef NASA_7_POLYNOMIAL_H
#define NASA_7_POLYNOMIAL_H

#include <iostream>

namespace Mutation {
    namespace Thermodynamics {


/**
 * Base class for NASA thermodynamic polynomial types.  Impliments polynomial
 * and temperature range storage.
 *
 * @todo Nasa9Polynomial does not actually derive from this so I should figure
 * out how make this nicer.
 */
template <int NR, int NC>
class NasaPolynomial
{
public:
    
    NasaPolynomial()
        : m_nr(NR), m_nc(NC)
    {
        for (int j = 0; j < NR; j++) {
            mp_tbounds[j] = 0.0;
            for (int i = 0; i < NC; i++)
                mp_coefficients[j][i] = 0.0;
        }
        mp_tbounds[NR] = 0.0;
    }
    
    NasaPolynomial(const NasaPolynomial& poly)
        : m_nr(NR), m_nc(NC)
    {
        for (int j = 0; j < NR; j++) {
            mp_tbounds[j] = poly.mp_tbounds[j];
            for (int i = 0; i < NC; i++)
                mp_coefficients[j][i] = poly.mp_coefficients[j][i];
        }
        mp_tbounds[NR] = poly.mp_tbounds[NR];
    }
    
    virtual ~NasaPolynomial() {};
    
    double minT() const {
        return mp_tbounds[0];
    }
    
    double maxT() const {
        return mp_tbounds[NR];
    }
    
    int tRange(double T) const
    {
        for (int i = 1; i < NR; ++i)
            if (T < mp_tbounds[i]) return i-1;
        return NR-1;
    }
    
    static int nCoefficients() { return NC; }
    
protected:
    
    const int m_nr;
    const int m_nc;
    
    double mp_coefficients[NR][NC];
    double mp_tbounds[NR+1];

};

/**
 * Represents a set of NASA 7-coefficient polynomials for computing
 * thermodynamic properties of pure substances.  The object consists of two 
 * polynomials that are defined for two separate but connected temperature 
 * ranges.
 */
class Nasa7Polynomial : public NasaPolynomial<2,7>
{
public:
    
    /**
     * Default constructor. Initializes all data members to 0.
     */
    Nasa7Polynomial() { };
    
    /**
     * Copy constructor.
     */
    Nasa7Polynomial(const Nasa7Polynomial &to_copy);
    
    /**
     * Assignment operator.
     */
    Nasa7Polynomial &operator=(const Nasa7Polynomial &to_copy);
    
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
     * Computes dimensionless specific heat, enthalpy, entropy, and gibbs free
     * energy - Cp/Ru, H/Ru/T, S/Ru, and G/Ru/T.
     * @see computeParams()
     */
    void thermo(const double *const p_params, double &cp, double &h, double &s,
                double &g) const;
    
    /**
     * Enumerates functions that require a parameters list.
     * @see computeParams()
     */
    enum ThermoFunction
    {
        CP,
        ENTHALPY,
        ENTROPY,
        GIBBS
    };
    
    /**
     * Computes the constant parameters as functions of temperature that must
     * be passed to the individual thermodynamic functions. The params array is
     * filled with the parameters needed by the function specified by func.  The
     * dimension of params should be at least as long as the following depending
     * on the value of func:
     *
     *  CP       - params[4]    \n
     *  ENTHALPY - params[6]    \n
     *  ENTROPY  - params[5]    \n
     *  GIBBS    - params[7]    
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

    friend std::istream &operator>>(std::istream &in, Nasa7Polynomial &n7);

};

/**
 * Insertion operator for Nasa7Polynomial object.  Initializes a Nasa7Polynomial
 * object from an input stream formatted in the manner specified by the NASA
 * database format for the 7 coefficient polynomial expressions for
 * thermodynamic properties.
 *
 * @see Nasa7Polynomial
 */
std::istream &operator>>(std::istream &in, Nasa7Polynomial &n7);

} // namespace Thermodynamics
} // namespace Mutation

#endif // NASA_7_POLYNOMIAL_H
