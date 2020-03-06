/**
 * @file bprime.cpp
 *
 * @brief Computes the solution of the surface mass balance equation and prints
 * a table to the console. @see @ref bprime.
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

#include "mutation++.h"
#include <iostream>
#include <Eigen/Dense>

#ifdef _GNU_SOURCE
#include <fenv.h>
#endif

using namespace std;
using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

/**
 * @page bprime B' Solver (bprime)
 * __Usage__:
 *
 *    bprime -T \f$T_1:\Delta T:T_2\f$ -p
 *    \f$p\f$ -b \f$B'_g\f$ -m mixture -bl BL -py Pyrolysis
 *
 * This program generates a so-called "B-prime" table for a given temperature
 * range and stepsize in K, a fixed pressure in Pa, a value of \f$B'_g\f$
 * (pyrolysis mass flux, nondimensionalized by the boundary layer edge mass
 * flux), a given mixture name.  Currently, this program assumes the char
 * composition to be solid graphite (which must be included in the mixture file).
 * The `BL` and `Pyrolysis` arguments are the names of the
 * [elemental compositions](@ref compositions) in the mixture file which
 * represent the boundary layer edge and pyrolysis gases respectively. The
 * produced table provides values of \f$B'_c\f$, the wall enthalpy in MJ/kg, and
 * the species mole fractions at the wall versus temperature.
 */
// Simply stores the command line options
typedef struct {
    double T1;
    double T2;
    double dT;

    double P1;
    double Bg;

    std::string mixture;
    std::string boundary_layer_comp;
    std::string pyrolysis_composition;

    bool pyrolysis_exist = false;
} Options;

// Checks if an option is present
bool optionExists(int argc, char** argv, const std::string& option)
{
    return (std::find(argv, argv+argc, option) != argv+argc);
}

// Returns the value associated with a particular option
std::string getOption(int argc, char** argv, const std::string& option)
{
    std::string value;
    char** ptr = std::find(argv, argv+argc, option);

    if (ptr == argv+argc || ptr+1 == argv+argc)
        value = "";
    else
        value = *(ptr+1);

    return value;
}


// Prints the program's usage information and exits.
void printHelpMessage(const char* const name)
{
    std::string tab("    ");

    cout.setf(std::ios::left, std::ios::adjustfield);

    cout << endl;
    cout << "Usage: " << name << " [OPTIONS] mixture" << endl;
    cout << "Compute the non-dimensional mass blowing rate for mixture "
         << "over a set of temperatures and pressure using the Mutation++ library." << endl;
    cout << endl;
    cout << tab << "-h, --help          prints this help message" << endl;
    cout << tab << "-T                  temperature range in K \"T1:dT:T2\" or simply T (default = 300:100:5000 K)" << endl;
    cout << tab << "-P                  pressure in Pa P (default = 1 atm)" << endl;
    cout << tab << "-b                  pyrolysis non-dimensional mass blowing rate (default = 0)" << endl;
    cout << tab << "-m                  mixture name" << endl;
    cout << tab << "-bl                 boundary layer edge composition name" << endl;
    cout << tab << "-py                 pyrolysis composition name (default = null)" << endl;


    cout << endl;
    cout << "Example:" << endl;
    cout << tab << name << " -T 300:100:5000 -P 101325 -b 10 -m carbonPhenol -bl BLedge -py Gas" << endl;
    cout << endl;
    cout << "Mixture file:" << endl;
    cout << tab << "carbonPhenol - corresponds to the name of the mixture" << endl;
    cout << tab << "BLedge - corresponds to the boundary layer edge elemental composition" << endl;
    cout << tab << "Gas - corresponds to the pyrolysis elemental gas composition" << endl;
    cout << endl;

    exit(0);
}

// Parses a temperature or pressure range
bool parseRange(const std::string& range, double& x1, double& x2, double& dx)
{
    std::vector<std::string> tokens;
    String::tokenize(range, tokens, ":");

    if (!String::isNumeric(tokens))
        return false;

    switch (tokens.size()) {
        case 1:
            x1 = atof(tokens[0].c_str());
            x2 = x1;
            dx = 1.0;
            break;
        case 3:
            x1 = atof(tokens[0].c_str());
            x2 = atof(tokens[2].c_str());
            dx = atof(tokens[1].c_str());
            break;
        default:
            return false;
    }

    if (dx == 0.0) {
        x2 = x1;
        dx = 1.0;
    }

    return true;
}


bool parseRange(const std::string& range, double& x1)
{
    std::vector<std::string> tokens;
    String::tokenize(range, tokens, ":");

    if (!String::isNumeric(tokens))
        return false;

    x1 = atof(tokens[0].c_str());


    return true;
}

// Parse the command line options to determine what the user wants to do
Options parseOptions(int argc, char** argv)
{
    Options opts;

    // Print the help message and exit if desired
    if(argc<2)
        printHelpMessage(argv[0]);

    if (optionExists(argc, argv, "-h") || optionExists(argc, argv, "--help"))
        printHelpMessage(argv[0]);

    // Get the temperature range
    if (optionExists(argc, argv, "-T")) {
        if (!parseRange(
                getOption(argc, argv, "-T"), opts.T1, opts.T2, opts.dT)) {
            cout << "Bad format for temperature range!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.T1 = 300.0;
        opts.T2 = 15000.0;
        opts.dT = 100.0;
    }

    // Get the pressure range
    if (optionExists(argc, argv, "-P")) {
        if (!parseRange(
                getOption(argc, argv, "-P"), opts.P1)) {
            cout << "Bad format for pressure !" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.P1 = ONEATM;
    }

    // Get the pressure range
    if (optionExists(argc, argv, "-b")) {
        if (!parseRange(
                getOption(argc, argv, "-b"), opts.Bg)) {
            cout << "Bad format for pressure !" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.Bg = 0;
    }

    if (optionExists(argc, argv, "-m")) {
        opts.mixture = getOption(argc, argv, "-m");
    } else {
        printHelpMessage(argv[0]);
    }


    if (optionExists(argc, argv, "-bl")) {
        opts.boundary_layer_comp = getOption(argc, argv, "-bl");
    } else {
        printHelpMessage(argv[0]);
    }

    if (optionExists(argc, argv, "-py")) {
        opts.pyrolysis_composition = getOption(argc, argv, "-py");
        opts.pyrolysis_exist= true;
    }

    return opts;
}


int main(int argc, char* argv[])
{
#ifdef _GNU_SOURCE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    Options opts = parseOptions(argc, argv);

    Mixture mix(opts.mixture);
    
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    
    std::vector<double> Yke (ne,0);
    std::vector<double> Ykg (ne,0);
    std::vector<double> Xw (ns,0);

    // Run conditions
    double T1 = opts.T1;
    double T2 = opts.T2;
    double dt = opts.dT;
    double P  = opts.P1;
    double Bg = opts.Bg;
    double Bc, hw;
    
    mix.getComposition(opts.boundary_layer_comp, Yke.data(), Composition::MASS);

    if(opts.pyrolysis_exist)
        mix.getComposition(opts.pyrolysis_composition, Ykg.data(), Composition::MASS);

    cout << setw(10) << "\"Tw[K]\""
         << setw(15) << "\"B'c\""
         << setw(15) << "\"hw[MJ/kg]\"";
    for (int i = 0; i < ns; ++i)
        cout << setw(25) << "\"" + mix.speciesName(i) + "\"";
    cout << endl;
    
    for (double T = T1; T < T2 + 1.0e-6; T += dt) {
        mix.surfaceMassBalance(Yke.data(), Ykg.data(), T, P, Bg, Bc, hw, Xw.data());
        cout << setw(10) << T << setw(15) << Bc << setw(15) << hw / 1.0e6;
        for (int i = 0; i < ns; ++i)
            cout << setw(25) << Xw[i];
        cout << endl;
    }

}
