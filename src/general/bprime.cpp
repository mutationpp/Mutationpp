/**
 * @file bprime.cpp
 *
 * @brief Computes the solution of the surface mass balance equation and prints
 * a table to the console. @see @ref bprime.
 */

/*
 * Copyright 2014-2017 von Karman Institute for Fluid Dynamics (VKI)
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

#include "mutation++.h"
#include <iostream>

#ifdef _GNU_SOURCE
#include <fenv.h>
#endif

using namespace std;
using namespace Mutation;
using namespace Mutation::Thermodynamics;

/**
 * @page bprime B' Solver (bprime)
 * __Usage__:
 *
 *    bprime \f$T_1\f$ \f$T_2\f$ \f$\Delta T\f$ \f$p\f$ \f$B'_g\f$ mixture BL Pyrolysis
 *
 * This program generates a so-called "B-prime" table for a given temperature
 * range and stepsize in K, a fixed pressure in atm, a value of \f$B'_g\f$
 * (pyrolysis mass flux, nondimensionalized by the boundary layer edge mass
 * flux), a given mixture name.  Currently, this program assumes the char
 * composition to be solid graphite (which must be included in the mixture file).
 * The `BL` and `Pyrolysis` arguments are the names of the
 * [elemental compositions](@ref compositions) in the mixture file which
 * represent the boundary layer edge and pyrolysis gases respectively. The
 * produced table provides values of \f$B'_c\f$, the wall enthalpy in MJ/kg, and
 * the species mole fractions at the wall versus temperature.
 */

int main(int argc, char* argv[])
{
#ifdef _GNU_SOURCE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    if (argc != 9) {
        cout << "Usage: bprime T1 T2 dT P bg mixture BL_comp Pyro_comp" << endl;
        exit(1);
    }

    Mixture mix(argv[6]);
    
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    double* p_Yke = new double [ne];
    double* p_Ykg = new double [ne];
    double* p_Xw  = new double [ns];
    
    // Run conditions
    double T1 = atof(argv[1]);
    double T2 = atof(argv[2]);
    double dt = atof(argv[3]);
    double P  = atof(argv[4])*ONEATM;
    double Bg = atof(argv[5]);
    double Bc, hw;
    
    mix.getComposition(argv[7], p_Yke, Composition::MASS);
    mix.getComposition(argv[8], p_Ykg, Composition::MASS);
    
    cout << setw(10) << "\"Tw[K]\""
         << setw(15) << "\"B'c\""
         << setw(15) << "\"hw[MJ/kg]\"";
    for (int i = 0; i < ns; ++i)
        cout << setw(25) << "\"" + mix.speciesName(i) + "\"";
    cout << endl;
    
    for (double T = T1; T < T2 + 1.0e-6; T += dt) {
        mix.surfaceMassBalance(p_Yke, p_Ykg, T, P, Bg, Bc, hw, p_Xw);
        cout << setw(10) << T << setw(15) << Bc << setw(15) << hw / 1.0e6;
        for (int i = 0; i < ns; ++i)
            cout << setw(25) << p_Xw[i];
        cout << endl;
    }
    
    delete [] p_Yke;
    delete [] p_Ykg;
    delete [] p_Xw;
}
