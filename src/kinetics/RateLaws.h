/**
 * @file RateLaws.h
 *
 * @brief Declaration of various RateLaw classes.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef RATELAW_H
#define RATELAW_H

#include <vector>
#include <cmath>
#include <cstdlib>

#include "Utilities.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Abstract base class for all rate laws which allows owners such as class 
 * Reaction to store any rate law polymorphically.
 */
class RateLaw
{
public:

    virtual ~RateLaw() { };
    virtual RateLaw* clone() const = 0;
};

/**
 * Arrhenius rate law \f$ k_f(T) = A T^\eta exp(-E_a / (R_u T)) \f$.
 */
class Arrhenius : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    Arrhenius(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    Arrhenius(const Arrhenius& to_copy)
        : m_lnA(to_copy.m_lnA), m_n(to_copy.m_n), m_temp(to_copy.m_temp)
    { }
    
    virtual ~Arrhenius() { };
    
    Arrhenius* clone() const {
        return new Arrhenius(*this);
    }
    
    inline double getLnRate(const double lnT, const double invT) const {
        return (m_lnA + m_n * lnT - m_temp * invT);
    }
    
    double A() const { 
        return std::exp(m_lnA);
    }
    
    double n() const {
        return m_n;
    }
    
    double T() const {
        return m_temp;
    }
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;    
    static std::vector<Mutation::Utilities::Units> sm_eunits;

    double m_lnA;
    double m_n;
    double m_temp;
};

class O2vO2 : public RateLaw {
public:
    O2vO2 ( const Mutation::Utilities::IO::XmlElement& node );
    ~O2vO2(){ }

    O2vO2( const O2vO2& copy) : m_vib_level(copy.m_vib_level){}
    O2vO2* clone() const { return new O2vO2(*this); }

    inline double getLnRate( const double lnT, const double invT ) const {

        const double l_lna = 6 * log(10.E0) + lnT - std::log( 1.8 ) - 12.E0 * std::log(10) - 122.E0 * std::pow( invT, 1.E0/3.E0 ) - std::log( 1.E0 - std::exp( - 2273.7E0 * invT )) ;
        const double l_b = 2.99E0 * std::sqrt( invT );

        return ( std::log(m_vib_level) + l_lna + l_b * ( m_vib_level - 1 ) );
    }

    double A() const { 
        return 0.E0;
    }
    
    double n() const {
        return 0.E0;
    }
    
    double T() const {
        return 0.E0;
    }

private:
    int m_vib_level;

};

class O2vO2w : public RateLaw {
public:
    O2vO2w ( const Mutation::Utilities::IO::XmlElement& node );
    ~O2vO2w(){ }
 
    O2vO2w( const O2vO2w& copy) : m_vib_level_v(copy.m_vib_level_v), 
                                  m_vib_level_w(copy.m_vib_level_w){}
    O2vO2w* clone() const { return new O2vO2w(*this); }

    inline double getLnRate( const double lnT, const double invT ) const {

	return ( std::log(2.8) - 24 * std::log(10) + std::log( m_vib_level_v ) + std::log(m_vib_level_w) + 1.5E0 * lnT + 2.4E0 * std::sqrt( invT ) * ( (m_vib_level_v - 1) - ( m_vib_level_w - 1 ) ) );
    }

    double A() const { 
	return 0.E0;
    }
    
    double n() const {
	return 0.E0;
    }
    
    double T() const {
	return 0.E0;
    }

private:
    int m_vib_level_v;
    int m_vib_level_w;

};

    } // namespace Kinetics
} // namespace Mutation

#endif // RATELAW_H
