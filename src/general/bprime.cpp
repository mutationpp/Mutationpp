/**
 * @file bprime.cpp
 *
 * @brief Computes the solution of the surface mass balance equation and prints
 * a table to the console. @see @ref bprime.
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


    if (argc != 10) {
        cout << "Usage: bprime T1 T2 dT P bg mixture BLedge Pyro Char" << endl;
        exit(1);
    }


   /* string mixture_name = argv[6];
    string wall_name = "_wall";
    string total_mixture_name = mixture_name+wall_name;
            cout<<total_mixture_name<<endl;*/

    Mixture mix(argv[6]);


    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    double* p_Ykc = new double [ne];
    double* p_Yke = new double [ne];
    double* p_Ykg = new double [ne];
    double* p_Xw  = new double [ns];
    double* p_Xkw  = new double [ne];




    
    // Run conditions
    double T1 = atof(argv[1]);
    double T2 = atof(argv[2]);
    double dt = atof(argv[3]);
    double P  = atof(argv[4])*ONEATM;
    double Bg = atof(argv[5]);
    double Bc, hw;

    
    mix.getComposition(argv[7], p_Yke, Composition::MASS);
    mix.getComposition(argv[8] , p_Ykg, Composition::MASS);
    mix.getComposition(argv[9] , p_Ykc, Composition::MASS);



 #ifdef Debug

    for(int i=0; i<ne; ++i){
        cout<<"p_Yke("<<mix.elementName(i)<<"): " <<p_Yke[i]<<endl;
        }
    for(int i=0; i<ne; ++i){
        cout<<"p_Ykg("<<mix.elementName(i)<<"): " <<p_Ykg[i]<<endl;
        }
    for(int i=0; i<ne; ++i){
        cout<<"p_Ykc("<<mix.elementName(i)<<"): " <<p_Ykc[i]<<endl;
        }
    cout<<"Bg: "<<Bg<<endl;

 #endif
   cout << setw(10) << "Variables =" <<setw(10) << "\"Tw[K]\""
         << setw(15) << "\"B'c\""
         << setw(15) << "\"hw[MJ/kg]\"";
    for (int i = 0; i < ns; ++i){
        cout << setw(15) << "\"" << mix.speciesName(i) << "\"";}
    for (int i = 0; i < ne; ++i){
        cout << setw(15) << "\"element_" << mix.elementName(i) << "\"";}
    cout << endl;

   /* cout<< "Variables =" <<setw(10)<<"\"Tw[K]\"" <<setw(15);
    for (int i = 0; i<ne;++i){
    cout<<"\""<<mix.elementName(i)<<"\""<<setw(15);
    }cout<<endl;*/



    for (double T = T1; T < T2 + 1.0e-6; T += dt) {
        for (int i = 0; i<ne; ++i){
            p_Xkw[i] = 0;
        }
        mix.surfaceMassBalance(p_Ykc, p_Yke, p_Ykg, T, P, Bg, Bc, hw, p_Xkw, p_Xw);

        cout << setw(10) << T << setw(15) << Bc << setw(15) << hw / 1.0e6;
        for (int i = 0; i < ns; ++i){
            cout << setw(15) << p_Xw[i];}
        for (int i = 0; i < ne; ++i){
            cout << setw(15) << p_Xkw[i];}

        cout << endl;
    }
    delete [] p_Xkw;
    delete [] p_Yke;
    delete [] p_Ykg;
    delete [] p_Xw;
    delete [] p_Ykc;
}
