#ifndef GENERAL_CONSTANTS_H
#define GENERAL_CONSTANTS_H

#include <cmath>

// Pi
const double PI    = 4.0 * std::atan(1.0);
const double TWOPI = 2.0 * PI;

// Universal constants (all units in SI)
const double NA     = 6.0221415E23;    // Avagadro's number (molecule/mol)
const double KB     = 1.3806503E-23;   // Boltzmann's constant (J/molecule-K)
const double RU     = NA * KB;         // Universal Gas constant (J/mole-K)
const double HP     = 6.626068E-34;    // Planck's constant (J-s)
const double MU0    = PI * 4.0E-7;     // Magnetic constant (H/m)
const double C0     = 299792458.0;     // Speed of light in vacuum (m/s)
const double EPS0   = 1.0/(MU0*C0*C0); // Vacuum permittivity (F/m)
const double QE     = 1.602176565E-19; // Elementary positive charge (C)
const double ONEATM = 101325.0;        // 1 atm in Pa
const double SQRT2  = std::sqrt(2.0);  // square root of 2

#endif // GENERAL_CONSTANTS_H
