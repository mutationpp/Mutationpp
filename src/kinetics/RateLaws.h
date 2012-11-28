#ifndef RATELAW_H
#define RATELAW_H

#include <vector>
#include <cmath>
#include <cstdlib>

#include "Units.h"
#include "XMLite.h"

using namespace utils;

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

    static void setUnits(XmlElement& node);
    
    Arrhenius(XmlElement& node, const int order);
    
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

    static std::vector<Units> sm_aunits;    
    static std::vector<Units> sm_eunits;

    double m_lnA;
    double m_n;
    double m_temp;
};


#endif // RATELAW_HPP
