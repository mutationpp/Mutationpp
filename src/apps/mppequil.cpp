/**
 * @file mppequil.cpp
 *
 * @brief Utility which provides equilibrium properties for a given mixture.
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
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

#include <Eigen/Dense>
using namespace Eigen;

#ifdef _GNU_SOURCE
#include <fenv.h>
#endif

using std::cout;
using std::endl;
using std::setw;

using namespace Mutation;
using namespace Mutation::Utilities;

#define COLUMN_WIDTH 14

// Defines a value that can be printed
struct OutputQuantity {
    std::string name;
    std::string units;
    std::string description;

    OutputQuantity(
        const std::string& n, const std::string& u, const std::string& d)
        : name(n), units(u), description(d)
    { }
};

// List of all mixture output quantities
#define NMIXTURE 45
OutputQuantity mixture_quantities[NMIXTURE] = {
    OutputQuantity("Th", "K", "heavy particle temperature"),
    OutputQuantity("P", "Pa", "pressure"),
    OutputQuantity("B", "T", "magnitude of magnetic field"),
    OutputQuantity("rho", "kg/m^3", "density"),
    OutputQuantity("nd", "1/m^3", "number density"),
    OutputQuantity("Mw", "kg/mol", "molecular weight"),
    OutputQuantity("Cp_eq", "J/mol-K", "equilibrium specific heat at constant pressure"),
    OutputQuantity("H", "J/mol", "mixture enthalpy"),
    OutputQuantity("S", "J/mol-K", "entropy"),
    OutputQuantity("Cp_eq", "J/kg-K", "equilibrium specific heat at constant pressure"),
    OutputQuantity("H", "J/kg", "mixture enthalpy"),
    OutputQuantity("H-H0", "J/kg", "mixture enthalpy minus the enthalpy at 0K"),
    OutputQuantity("S", "J/kg-K", "entropy"),
    OutputQuantity("Cv_eq", "J/kg-K", "equilibrium specific heat at constant volume"),
    OutputQuantity("Cp", "J/mol-K", "frozen specific heat at constant pressure"),
    OutputQuantity("Cv", "J/mol-K", "frozen specific heat at constant volume"),
    OutputQuantity("Cp", "J/kg-K", "frozen specific heat at constant pressure"),
    OutputQuantity("Cv", "J/kg-K", "frozen specific heat at constant volume"),
    OutputQuantity("gam_eq", "", "equilibrium ratio of specific heats"),
    OutputQuantity("gamma", "", "frozen ratio of specific heat"),
    OutputQuantity("Ht", "J/mol", "translational enthalpy"),
    OutputQuantity("Hr", "J/mol", "rotational enthalpy"),
    OutputQuantity("Hv", "J/mol", "vibrational enthalpy"),
    OutputQuantity("Hel", "J/mol", "electronic enthalpy"),
    OutputQuantity("Hf", "J/mol", "formation enthalpy"),
    OutputQuantity("Ht", "J/kg", "translational enthalpy"),
    OutputQuantity("Hr", "J/kg", "rotational enthalpy"),
    OutputQuantity("Hv", "J/kg", "vibrational enthalpy"),
    OutputQuantity("Hel", "J/kg", "electronic enthalpy"),
    OutputQuantity("Hf", "J/kg", "formation enthalpy"),
    OutputQuantity("e", "J/mol", "mixture energy"),
    OutputQuantity("e", "J/kg", "mixture energy"),
    OutputQuantity("mu", "Pa-s", "dynamic viscosity"),
    OutputQuantity("lambda", "W/m-K", "mixture equilibrium thermal conductivity"),
    OutputQuantity("lam_reac", "W/m-K", "reactive thermal conductivity"),
    OutputQuantity("lam_bb", "W/m-K", "Butler-Brokaw reactive thermal conductivity"),
    OutputQuantity("lam_soret", "W/m-K", "Soret thermal conductivity"),
    OutputQuantity("lam_int", "W/m-K", "internal energy thermal conductivity"),
    OutputQuantity("lam_h", "W/m-K", "heavy particle translational thermal conductivity"),
    OutputQuantity("lam_e", "W/m-K", "electron translational thermal conductivity"),
    OutputQuantity("sigma", "S/m", "electric conductivity (B=0)"),
//    OutputQuantity("sigma_para", "S/m", "electronic conductivity parallel to magnetic field"),
//    OutputQuantity("sigma_perp", "S/m", "electronic conductivity perpendicular to magnetic field"),
//    OutputQuantity("sigma_tran", "S/m", "electronic conductivity transverse to magnetic field"),
    OutputQuantity("a_f", "m/s", "frozen speed of sound"),
    OutputQuantity("a_eq", "m/s", "equilibrium speed of sound"),
    OutputQuantity("Eam", "V/K", "ambipolar electric field (SM Ramshaw)"),
    OutputQuantity("drho/dP", "kg/J", "equilibrium density derivative w.r.t pressure")
//    OutputQuantity("l", "m", "mean free path"),
//    OutputQuantity("le", "m", "mean free path of electrons"),
//    OutputQuantity("Vh", "m/s", "average heavy particle thermal speed"),
//    OutputQuantity("Ve", "m/s", "electron thermal speed"),
//    OutputQuantity("tau_eh", "1/s", "electron-heavy collision frequency"),
//    OutputQuantity("tau_a", "1/s", "average heavy particle collision frequency")
};

// List of all species output quantities
#define NSPECIES 19
OutputQuantity species_quantities[NSPECIES] = {
    OutputQuantity("X", "", "mole fractions"),
    OutputQuantity("dX/dT", "1/K", "partial of mole fraction w.r.t. temperature"),
    OutputQuantity("Y", "", "mass fractions"),
    OutputQuantity("rho", "kg/m^3", "mass densities"),
    OutputQuantity("conc", "mol/m^3", "molar concentrations"),
    OutputQuantity("Cp", "J/mol-K", "specific heats at constant pressure"),
    OutputQuantity("H", "J/mol", "enthalpies"),
    OutputQuantity("S", "J/mol-K", "entropies"),
    OutputQuantity("G", "J/mol", "Gibbs free energies"),
    OutputQuantity("Cp", "J/kg-K", "specific heats at constant pressure"),
    OutputQuantity("H", "J/kg", "enthalpies"),
    OutputQuantity("S", "J/kg-K", "entropies"),
    OutputQuantity("G", "J/kg", "Gibbs free energies"),
    OutputQuantity("J", "kg/m^2-s", "Species diffusion fluxes (SM Ramshaw)"),
    OutputQuantity("omega", "kg/m^3-s", "production rates due to reactions"),
    OutputQuantity("Omega11", "m^2", "(1,1) pure species collision integrals"),
    OutputQuantity("Omega22", "m^2", "(2,2) pure species collision integrals"),
    OutputQuantity("Chi^h", "", "heavy thermal diffusion ratios"),
    OutputQuantity("Dm", "m^2/s", "mixture averaged diffusion coefficients")
};

// List of reaction output quantities
#define NREACTION 2
OutputQuantity reaction_quantities[NREACTION] = {
    OutputQuantity("kf", "mol,m,s,K", "forward reaction rate coefficients"),
    OutputQuantity("kb", "mol,m,s,K", "backward reaction rate coefficients")
};

// List of other output quantities
#define NOTHER 10
OutputQuantity other_quantities[NOTHER] = {
    OutputQuantity("Dij", "m^2/s", "multicomponent diffusion coefficients"),
    OutputQuantity("pi_i", "", "element potentials"),
    OutputQuantity("N_p", "mol", "phase moles"),
    OutputQuantity("iters", "", "number of continuation step iterations"),
    OutputQuantity("newts", "", "total number of newton iterations"),
    OutputQuantity("Fp_k", "kg/m-Pa-s", "elemental diffusion fluxes per pressure gradient"),
    OutputQuantity("Ft_k", "kg/m-K-s", "elemental diffusion fluxes per temperature gradient"),
    OutputQuantity("Fz_k", "kg/m-s", "elemental diffusion fluxes per element mole fraction gradient"),
    OutputQuantity("sigmaB", "S/m", "anisotropic electric conductivity"),
    OutputQuantity("lamB_e", "W/m-K", "anisotropic electron thermal conductivity")
};

// Simply stores the command line options
typedef struct {
    double T1;
    double T2;
    double dT;

    double P1;
    double P2;
    double dP;
    double B;

    std::vector<int> mixture_indices;
    std::vector<int> species_indices;
    std::vector<int> reaction_indices;
    std::vector<int> other_indices;

    bool verbose;
    bool header;

    Mutation::MixtureOptions* p_mixture_opts;

    bool use_scientific;
    int  precision;
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
    cout << "Compute equilibrium properties for mixture over a set of "
         << "temperatures and pressures using the Mutation++ library." << endl;
    cout << endl;
    cout << tab << "-h, --help          prints this help message" << endl;
    cout << tab << "    --no-header     no table header will be printed" << endl;
    cout << tab << "-T                  temperature range in K \"T1:dT:T2\" or simply T (default = 300:100:15000 K)" << endl;
    cout << tab << "-P                  pressure range in Pa \"P1:dP:P2\" or simply P (default = 1 atm)" << endl;
    cout << tab << "-B                  magnitude of the magnetic field in teslas (default = 0 T)" << endl;
    cout << tab << "-m                  list of mixture values to output (see below)" << endl;
    cout << tab << "-s                  list of species values to output (see below)" << endl;
    cout << tab << "-r                  list of reaction values to output (see below)" << endl;
    cout << tab << "-o                  list of other values to output (see below)" << endl;
    cout << tab << "    --species-list  instead of mixture name, use this to list species in mixture" << endl;
    cout << tab << "    --elem-x        set elemental mole fractions (ex: N:0.8,O:0.2)" << endl;
    cout << tab << "    --elem-comp     set elemental composition with a name from the mixture file" << endl;
    cout << tab << "    --thermo-db     overrides thermodynamic database type (NASA-7, NASA-9, RRHO)" << endl;
    cout << tab << "    --scientific    outputs in scientific format with given precision" << endl;
    //cout << tab << "-c             element fractions (ie: \"N:0.79,O:0.21\")" << endl;
    cout << endl;
    cout << "Mixture values (example format: \"1-3,7,9-11\"):" << endl;

    for (int i = 0; i < NMIXTURE; ++i)
        cout << tab << setw(2) << i << ": " << setw(10)
             << mixture_quantities[i].name << setw(12)
             << (mixture_quantities[i].units == "" ? "[-]" :
                "[" + mixture_quantities[i].units + "]")
             << mixture_quantities[i].description << endl;

    cout << endl;
    cout << "Species values (same format as mixture values):" << endl;

    for (int i = 0; i < NSPECIES; ++i)
        cout << tab << setw(2) << i << ": " << setw(10)
             << species_quantities[i].name << setw(12)
             << (species_quantities[i].units == "" ? "[-]" :
                "[" + species_quantities[i].units + "]")
             << species_quantities[i].description << endl;

    cout << endl;
    cout << "Reaction values (same format as mixture values):" << endl;

    for (int i = 0; i < NREACTION; ++i)
        cout << tab << setw(2) << i << ": " << setw(10)
             << reaction_quantities[i].name << setw(12)
             << (reaction_quantities[i].units == "" ? "[-]" :
                "[" + reaction_quantities[i].units + "]")
             << reaction_quantities[i].description << endl;

    cout << endl;
    cout << "Other values (same format as mixture values):" << endl;

    for (int i = 0; i < NOTHER; ++i)
        cout << tab << setw(2) << i << ": " << setw(10)
             << other_quantities[i].name << setw(12)
             << (other_quantities[i].units == "" ? "[-]" :
                "[" + other_quantities[i].units + "]")
             << other_quantities[i].description << endl;

    cout << endl;
    cout << "Example:" << endl;
    cout << tab << name << " -T 300:100:15000 -P 101325 -m 1-3,8 air11" << endl;
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

// Parses an index list of the form 1-3,7,9 for example
bool parseIndices(const std::string& list, std::vector<int>& indices, int max)
{
    indices.clear();

    std::vector<std::string> ranges;
    std::vector<std::string> bounds;
    String::tokenize(list, ranges, ",");

    std::vector<std::string>::const_iterator iter = ranges.begin();
    for ( ; iter != ranges.end(); ++iter) {
        bounds.clear();
        String::tokenize(*iter, bounds, "-");

        if (!String::isNumeric(bounds))
            return false;

        switch (bounds.size()) {
            case 1: {
                int i = atoi(bounds[0].c_str());
                if (i < 0 || i > max)
                    return false;
                indices.push_back(atoi(bounds[0].c_str()));
                break;
            }
            case 2: {
                int i1 = atoi(bounds[0].c_str());
                int i2 = atoi(bounds[1].c_str());

                if (i1 >= i2 || i1 < 0 || i1 > max || i2 > max)
                    return false;

                for (int i = i1; i <= i2; ++i)
                    indices.push_back(i);

                break;
            }
            default:
                return false;
        }
    }

    return true;
}

// Parse the command line options to determine what the user wants to do
Options parseOptions(int argc, char** argv)
{
    Options opts;

    // Print the help message and exit if desired
    if (optionExists(argc, argv, "-h") || optionExists(argc, argv, "--help"))
        printHelpMessage(argv[0]);

    // Control verbosity
    opts.verbose =
        optionExists(argc, argv, "-v") || optionExists(argc, argv, "--verbose");

    // Check if the header should be printed or not
    opts.header = !optionExists(argc, argv, "--no-header");

    // The mixture name is given as the only argument (unless --species option
    // is present)
    if (optionExists(argc, argv, "--species-list")) {
        opts.p_mixture_opts = new Mutation::MixtureOptions();
        opts.p_mixture_opts->setSpeciesDescriptor(getOption(argc, argv, "--species-list"));
    } else {
        opts.p_mixture_opts = new Mutation::MixtureOptions(argv[argc-1]);
    }

    // Must use the equilibirum state model
    opts.p_mixture_opts->setStateModel("Equil");

    // Elemental mole fractions
    bool comp_set = false;
    if (optionExists(argc, argv, "--elem-x")) {
        comp_set = true;
        opts.p_mixture_opts->addComposition(
            Thermodynamics::Composition(getOption(argc, argv, "--elem-x").c_str()),
            true
        );
    }

    // Elemental composition
    if (optionExists(argc, argv, "--elem-comp")) {
        if (comp_set) {
            cout << "Only one method of setting element composition can be used!" << endl;
            exit(1);
        }
        std::string comp = getOption(argc, argv, "--elem-comp");
        if (!opts.p_mixture_opts->setDefaultComposition(comp.c_str())) {
            cout << "Composition " << comp << " does not exist in the mixture!" << endl;
            exit(1);
        }
    }

    // Thermodynamic database
    if (optionExists(argc, argv, "--thermo-db")) {
        opts.p_mixture_opts->setThermodynamicDatabase(
            getOption(argc, argv, "--thermo-db"));
    }

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
            getOption(argc, argv, "-P"), opts.P1, opts.P2, opts.dP)) {
            cout << "Bad format for pressure range!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.P1 = ONEATM;
        opts.P2 = ONEATM;
        opts.dP = ONEATM;
    }

    // Get the magnetic field strength
    if (optionExists(argc, argv, "-B")) {
        opts.B = atof(getOption(argc, argv, "-B").c_str());
    } else {
        opts.B = 0.0;
    }

    // Get the mixture properties to print
    if (optionExists(argc, argv, "-m")) {
        if (!parseIndices(
            getOption(argc, argv, "-m"), opts.mixture_indices, NMIXTURE-1)) {
            cout << "Bad format for mixture value indices!" << endl;
            printHelpMessage(argv[0]);
        }
    }

    // Get the species properties to print
    if (optionExists(argc, argv, "-s")) {
        if (!parseIndices(
            getOption(argc, argv, "-s"), opts.species_indices, NSPECIES-1)) {
            cout << "Bad format for species value indices!" << endl;
            printHelpMessage(argv[0]);
        }
    }

    // Get the reaction properties to print
    if (optionExists(argc, argv, "-r")) {
        if (!parseIndices(
            getOption(argc, argv, "-r"), opts.reaction_indices, NREACTION-1)) {
            printHelpMessage(argv[0]);
        }
    }

    // Get the other properties to print
    if (optionExists(argc, argv, "-o")) {
        if (!parseIndices(
            getOption(argc, argv, "-o"), opts.other_indices, NOTHER-1)) {
            cout << "Bad format for other value indices!" << endl;
            printHelpMessage(argv[0]);
        }
    }

    // Check for output format
    if (optionExists(argc, argv, "--scientific")) {
        opts.use_scientific = true;
        opts.precision = atoi(getOption(argc, argv, "--scientific").c_str());
    } else {
        opts.use_scientific = false;
        opts.precision      = 0;
    }

    // If

    return opts;
}

// Write out the column headers
void writeHeader(
    const Options& opts, const Mutation::Mixture& mix,
    std::vector<int>& column_widths)
{
    std::string name;
    int width = (opts.use_scientific ? opts.precision + 8 : COLUMN_WIDTH);

    std::vector<int>::const_iterator iter = opts.mixture_indices.begin();
    for ( ; iter != opts.mixture_indices.end(); ++iter) {
        name = mixture_quantities[*iter].name +
            (mixture_quantities[*iter].units == "" ?
                "" : "[" + mixture_quantities[*iter].units + "]");
        column_widths.push_back(
            std::max(width, static_cast<int>(name.length())+2));
        if (opts.header)
            cout << setw(column_widths.back()) << name;
    }

    iter = opts.species_indices.begin();
    for ( ; iter != opts.species_indices.end(); ++iter) {
        for (int i = 0; i < mix.nSpecies(); ++i) {
            name = "\"" + species_quantities[*iter].name + "_" + mix.speciesName(i) +
                (species_quantities[*iter].units == "" ?
                    "" : "[" + species_quantities[*iter].units + "]") + "\"";
            column_widths.push_back(
                std::max(width, static_cast<int>(name.length())+2));
            if (opts.header)
                cout << setw(column_widths.back()) << name;
        }
    }

    iter = opts.reaction_indices.begin();
    for ( ; iter != opts.reaction_indices.end(); ++iter) {
        for (int i = 0; i < mix.nReactions(); ++i) {
            char buff [10];
            sprintf(buff, "%d", i);
            name = reaction_quantities[*iter].name + "_" + buff +
                (reaction_quantities[*iter].units == "" ?
                    "" : "[" + reaction_quantities[*iter].units + "]");
            column_widths.push_back(
                std::max(width, static_cast<int>(name.length())+2));
            if (opts.header)
                cout << setw(column_widths.back()) << name;
        }
    }

    iter = opts.other_indices.begin();
    for ( ; iter != opts.other_indices.end(); ++iter) {
        if (other_quantities[*iter].name == "Dij") {
            for (int i = 0; i < mix.nSpecies(); ++i) {
                for (int j = 0; j < mix.nSpecies(); ++j) {
                    name = "D_{" + mix.speciesName(i) + "," + mix.speciesName(j)
                        + "}";
                    column_widths.push_back(
                        std::max(width, static_cast<int>(name.length())+2));
                    if (opts.header)
                        cout << setw(column_widths.back()) << name;
                }
            }
        } else if (other_quantities[*iter].name == "pi_i") {
            for (int i = 0; i < mix.nElements(); ++i) {
                name = "pi_" + mix.elementName(i);
                column_widths.push_back(
                    std::max(width, static_cast<int>(name.length())+2));
                if (opts.header)
                    cout << setw(column_widths.back()) << name;
            }
        } else if (other_quantities[*iter].name == "N_p") {
            for (int i = 0; i < mix.nPhases(); ++i) {
                std::stringstream ss;
                ss << "N_p" << i; ss >> name;
                column_widths.push_back(
                    std::max(width, static_cast<int>(name.length())+2));
                if (opts.header)
                    cout << setw(column_widths.back()) << name;
            }
        } else if (other_quantities[*iter].name == "iters") {
            column_widths.push_back(width);
            if (opts.header)
                cout << setw(column_widths.back()) << "iters";
        } else if (other_quantities[*iter].name == "newts") {
            column_widths.push_back(width);
            if (opts.header)
                cout << setw(column_widths.back()) << "newts";
        } else if (other_quantities[*iter].name == "Fp_k") {
            for (int i = 0; i < mix.nElements(); ++i) {
                std::stringstream ss;
                ss << "Fp_" << mix.elementName(i) << "[kg/m-Pa-s]"; ss >> name;
                column_widths.push_back(
                    std::max(width, static_cast<int>(name.length())+2));
                if (opts.header)
                    cout << setw(column_widths.back()) << name;
            }
        } else if (other_quantities[*iter].name == "Ft_k") {
            for (int i = 0; i < mix.nElements(); ++i) {
                std::stringstream ss;
                ss << "Ft_" << mix.elementName(i) << "[kg/m-K-s]"; ss >> name;
                column_widths.push_back(
                    std::max(width, static_cast<int>(name.length())+2));
                if (opts.header)
                    cout << setw(column_widths.back()) << name;
            }
        } else if (other_quantities[*iter].name == "Fz_k") {
            for (int l = 0; l < mix.nElements(); ++l) {
                for (int k = 0; k < mix.nElements(); ++k) {
                    std::stringstream ss;
                    ss << "F(" << mix.elementName(l) << ")_" <<
                       mix.elementName(k) << "[kg/m-s]"; ss >> name;
                    column_widths.push_back(
                        std::max(width, static_cast<int>(name.length())+2));
                    if (opts.header)
                        cout << setw(column_widths.back()) << name;
                }
            }
        } else if (other_quantities[*iter].name == "sigmaB") {
            column_widths.push_back(std::max(width, 10));
            if (opts.header)
                cout << setw(column_widths.back()) << "sig_par";
            column_widths.push_back(std::max(width, 10));
            if (opts.header)
                cout << setw(column_widths.back()) << "sig_perp";
            column_widths.push_back(std::max(width, 10));
            if (opts.header)
                cout << setw(column_widths.back()) << "sig_tran";
        } else if (other_quantities[*iter].name == "lamB_e") {
            column_widths.push_back(std::max(width, 10));
            if (opts.header)
                cout << setw(column_widths.back()) << "lamB_par";
            column_widths.push_back(std::max(width, 10));
            if (opts.header)
                cout << setw(column_widths.back()) << "lamB_perp";
            column_widths.push_back(std::max(width, 10));
            if (opts.header)
                cout << setw(column_widths.back()) << "lamB_tran";
        }
    }

    if (opts.header)
        cout << endl;
}


/**
 * @page mppequil Mutation++ Equilibrium Properties (mppequil)
 *
 * __Usage:__
 *
 * mppequil [OPTIONS] mixture
 *
 * Compute equilibrium properties for mixture over a set of temperatures and
 * pressures using the Mutation++ library.  Use
 *
 *     mppequil -h
 *
 * for a full list of options.
 */


int main(int argc, char** argv)
{
#ifdef _GNU_SOURCE
    // Enable floating point exception handling
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    // Parse the command line options and load the mixture
    Options opts = parseOptions(argc, argv);
    Mutation::Mixture mix(*opts.p_mixture_opts);
    mix.setBField(opts.B);
    delete opts.p_mixture_opts;

    // Write out the column headers (and compute the column sizes)
    std::vector<int> column_widths;
    writeHeader(opts, mix, column_widths);

    // Now we can actually perform the computations
    std::vector<int>::const_iterator iter;
    double value;
    int cw;
    std::string name, units;

    double* species_values = new double [mix.nSpecies()];
    double* temp           = new double [
		std::max(mix.nSpecies(), mix.nElements()*mix.nElements())];
    double* sjac_exact     = new double [mix.nSpecies()*mix.nSpecies()];
    double* sjac_fd        = new double [mix.nSpecies()*mix.nSpecies()];
    double* temp2          = new double [mix.nSpecies()];
    double* reaction_values = new double [mix.nReactions()];

    /*ofstream jac_file("jac.dat");

    for (int i = 0; i < mix.nSpecies(); ++i)
        jac_file << setw(15) << mix.speciesName(i);
    jac_file << endl;*/

    if (opts.use_scientific) {
        cout.precision(opts.precision);
        cout << std::scientific;
    }

    for (double P = opts.P1; P <= opts.P2; P += opts.dP) {
        for (double T = opts.T1; T <= opts.T2; T += opts.dT) {
            // Compute the equilibrium composition
            mix.setState(&P, &T, 1);
            cw = 0;

            // Mixture properties
            iter = opts.mixture_indices.begin();
            for ( ; iter < opts.mixture_indices.end(); ++iter) {
                name  = mixture_quantities[*iter].name;
                units = mixture_quantities[*iter].units;

                if (name == "Th")
                    value = mix.T();
                else if (name == "P")
                    value = mix.P();
                else if (name == "B")
                    value = mix.getBField();
                else if (name == "rho")
                    value = mix.density();
                else if (name == "nd")
                    value = mix.numberDensity();
                else if (name == "Mw")
                    value = mix.mixtureMw();
                else if (name == "H") {
                    if (units == "J/mol")
                        value = mix.mixtureHMole();
                    else if (units == "J/kg")
                        value = mix.mixtureHMass();
                } else if (name == "H-H0") {
                    value = mix.mixtureHMinusH0Mass();
                } else if (name == "S") {
                    if (units == "J/mol-K")
                        value = mix.mixtureSMole();
                    else if (units == "J/kg-K")
                        value = mix.mixtureSMass();
                } else if (name == "Cp") {
                    if (units == "J/mol-K")
                        value = mix.mixtureFrozenCpMole();
                    else if (units == "J/kg-K")
                        value = mix.mixtureFrozenCpMass();
                } else if (name == "Cp_eq") {
                    if (units == "J/mol-K")
                        value = mix.mixtureEquilibriumCpMole();
                    else if (units == "J/kg-K")
                        value = mix.mixtureEquilibriumCpMass();
                } else if (name == "Cv") {
                    if (units == "J/mol-K")
                        value = mix.mixtureFrozenCvMole();
                    else if (units == "J/kg-K")
                        value = mix.mixtureFrozenCvMass();
                } else if (name == "Cv_eq") {
                    value = mix.mixtureEquilibriumCvMass();
                } else if (name == "gam_eq")
                    value = mix.mixtureEquilibriumGamma();
                else if (name == "gamma")
                    value = mix.mixtureFrozenGamma();
                else if (name == "mu")
                    value = mix.viscosity();
                else if (name == "lambda")
                    value = mix.equilibriumThermalConductivity();
                else if (name == "lam_reac")
                    value = mix.reactiveThermalConductivity();
                else if (name == "lam_bb")
                    value = mix.butlerBrokawThermalConductivity();
                else if (name == "lam_soret")
                    value = mix.soretThermalConductivity();
                else if (name == "lam_int")
                    value = mix.internalThermalConductivity(T);
                else if (name == "lam_h")
                    value = mix.heavyThermalConductivity();
                else if (name == "lam_e")
                    value = mix.electronThermalConductivity();
                else if (name == "sigma")
                    value = mix.electricConductivity();
//                else if (name == "sigma_para")
//                    value = mix.sigmaParallel();
//                else if (name == "sigma_perp")
//                    value = mix.sigmaPerpendicular();
//                else if (name == "sigma_trans")
//                    value = mix.sigmaTransverse();
                else if (name == "Ht") {
                    mix.speciesHOverRT(temp, species_values);
                    value = 0.0;
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        value += species_values[i] * RU * T * mix.X()[i];
                    if (units == "J/kg")
                        value /= mix.mixtureMw();
                } else if (name == "Hr") {
                    mix.speciesHOverRT(temp, NULL, species_values);
                    value = 0.0;
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        value += species_values[i] * RU * T * mix.X()[i];
                    if (units == "J/kg")
                        value /= mix.mixtureMw();
                } else if (name == "Hv") {
                    mix.speciesHOverRT(temp, NULL, NULL, species_values);
                    value = 0.0;
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        value += species_values[i] * RU * T * mix.X()[i];
                    if (units == "J/kg")
                        value /= mix.mixtureMw();
                } else if (name == "Hel") {
                    mix.speciesHOverRT(temp, NULL, NULL, NULL, species_values);
                    value = 0.0;
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        value += species_values[i] * RU * T * mix.X()[i];
                    if (units == "J/kg")
                        value /= mix.mixtureMw();
                } else if (name == "Hf") {
                    mix.speciesHOverRT(
                        temp, NULL, NULL, NULL, NULL, species_values);
                    value = 0.0;
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        value += species_values[i] * RU * T * mix.X()[i];
                    if (units == "J/kg")
                        value /= mix.mixtureMw();
                } else if (name == "e") {
                    if (units == "J/mol")
                        value = mix.mixtureEnergyMole();
                    else if (units == "J/kg")
                        value = mix.mixtureEnergyMass();
                } else if (name == "a_f")
                    value = mix.frozenSoundSpeed();
                else if (name == "a_eq")
                    value = mix.equilibriumSoundSpeed();
                else if (name == "Eam") {
                    mix.dXidT(temp);
                    mix.stefanMaxwell(temp, temp2, value);
                } else if (name == "drho/dP")
                    value = mix.dRhodP();
//                else if (name == "l")
//                    value = mix.meanFreePath();
//                else if (name == "le")
//                    value = mix.electronMeanFreePath();
//                else if (name == "Vh")
//                    value = mix.averageHeavyThermalSpeed();
//                else if (name == "Ve")
//                    value = mix.electronThermalSpeed();
//                else if (name == "tau_eh")
//                    value = mix.electronHeavyCollisionFreq();
//                else if (name == "tau_a")
//                    value = mix.averageHeavyCollisionFreq();

                cout << setw(column_widths[cw++]) << value;
            }

            // Species properties
            iter = opts.species_indices.begin();
            for ( ; iter < opts.species_indices.end(); ++iter) {
                name  = species_quantities[*iter].name;
                units = species_quantities[*iter].units;

                if (name == "X")
                    std::copy(mix.X(), mix.X()+mix.nSpecies(), species_values);
                else if (name == "dX/dT")
                    mix.dXidT(species_values);
                else if (name == "Y")
                    std::copy(mix.Y(), mix.Y()+mix.nSpecies(), species_values);
                else if (name == "rho") {
                    std::copy(mix.Y(), mix.Y()+mix.nSpecies(), species_values);
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        species_values[i] *= mix.density();
                } else if (name == "conc") {
                    double conc = mix.density() / mix.mixtureMw();
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        species_values[i] = mix.X()[i] * conc;
                } else if (name == "Cp") {
                    mix.speciesCpOverR(species_values);
                    if (units == "J/mol-K")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= RU;
                    else if (units == "J/kg-K")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= (RU / mix.speciesMw(i));
                } else if (name == "H") {
                    mix.speciesHOverRT(species_values);
                    if (units == "J/mol")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= (RU * T);
                    else if (units == "J/kg")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= (RU * T / mix.speciesMw(i));
                } else if (name == "S") {
                    mix.speciesSOverR(species_values);
                    if (units == "J/mol-K")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= RU;
                    else if (units == "J/kg-K")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= (RU / mix.speciesMw(i));
                } else if (name == "G") {
                    mix.speciesGOverRT(species_values);
                    if (units == "J/mol")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= (RU * T);
                    else if (units == "J/kg")
                        for (int i = 0; i < mix.nSpecies(); ++i)
                            species_values[i] *= (RU * T / mix.speciesMw(i));
                } else if (name == "J") {
                    double E, rho;
                    mix.dXidT(temp);
                    mix.stefanMaxwell(temp, species_values, E);
                    rho = mix.density();
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        species_values[i] *= rho * mix.Y()[i];
                } else if (name == "omega") {
                    mix.netProductionRates(species_values);
                } else if (name == "Omega11") {
                    int k = mix.hasElectrons() ? 1 : 0;
                    if (k == 1)
                        species_values[0] = mix.collisionDB().Q11ee();
                    Map<ArrayXd>(species_values+k, mix.nHeavy()) =
                        mix.collisionDB().Q11ii();
                } else if (name == "Omega22") {
                    int k = mix.hasElectrons() ? 1 : 0;
                    if (k == 1)
                        species_values[0] = mix.collisionDB().Q22ee();
                    Map<ArrayXd>(species_values+k, mix.nHeavy()) =
                        mix.collisionDB().Q22ii();
                } else if (name == "Chi^h") {
                    mix.heavyThermalDiffusionRatios(species_values);
                } else if (name == "Dm") {
                    mix.averageDiffusionCoeffs(species_values);
                }

                for (int i = 0; i < mix.nSpecies(); ++i)
                    cout << setw(column_widths[cw++]) << species_values[i];
            }

            // Reaction properties
            iter = opts.reaction_indices.begin();
            for ( ; iter < opts.reaction_indices.end(); ++iter) {
                name  = reaction_quantities[*iter].name;

                if (name == "kf")
                    mix.forwardRateCoefficients(reaction_values);
                else if (name == "kb")
                    mix.backwardRateCoefficients(reaction_values);

                for (int i = 0; i < mix.nReactions(); ++i)
                    cout << setw(column_widths[cw++]) << reaction_values[i];
            }

            // Other properties
            iter = opts.other_indices.begin();
            for ( ; iter < opts.other_indices.end(); ++iter) {
                name  = other_quantities[*iter].name;

                if (name == "Dij") {
                    const Eigen::MatrixXd& Dij = mix.diffusionMatrix();
                    for (int i = 0; i < mix.nSpecies(); ++i)
                        for (int j = 0; j < mix.nSpecies(); ++j)
                            cout << setw(column_widths[cw++]) << Dij(i,j);
                } else if (name == "pi_i") {
                    mix.elementPotentials(species_values);
                    for (int i = 0; i < mix.nElements(); ++i)
                        cout << setw(column_widths[cw++]) << species_values[i];
                } else if (name == "N_p") {
                    mix.phaseMoles(species_values);
                    for (int i = 0; i < mix.nPhases(); ++i)
                        cout << setw(column_widths[cw++]) << species_values[i];
                } else if (name == "iters") {
                    cout << setw(column_widths[cw++]) << mix.nEquilibriumSteps();
                } else if (name == "newts") {
                    cout << setw(column_widths[cw++]) << mix.nEquilibriumNewtons();
                } else if (name == "Fp_k") {
                    mix.equilDiffFluxFacsP(temp);
                    for (int i = 0; i < mix.nElements(); ++i)
                        cout << setw(column_widths[cw++]) << temp[i];
                } else if (name == "Ft_k") {
                    mix.equilDiffFluxFacsT(temp);
                    for (int i = 0; i < mix.nElements(); ++i)
                        cout << setw(column_widths[cw++]) << temp[i];
                } else if (name == "Fz_k") {
                    mix.equilDiffFluxFacsZ(temp);
                    for (int k = 0; k < mix.nElements()*mix.nElements(); ++k)
                        cout << setw(column_widths[cw++]) << temp[k];
                } else if (name == "sigmaB") {
                    Eigen::Vector3d sigma = mix.electricConductivityB();
                    for (int i = 0; i < 3; ++i)
                        cout << setw(column_widths[cw++]) << sigma(i);
                } else if (name == "lamB_e") {
                    Eigen::Vector3d lambda = mix.electronThermalConductivityB();
                    for (int i = 0; i < 3; ++i)
                        cout << setw(column_widths[cw++]) << lambda(i);
                }
            }

            cout << endl;

            // Compute the exact jacobian
            /*double conc = mix.density() / mix.mixtureMw();
            for (int i = 0; i < mix.nSpecies(); ++i)
                temp[i] = mix.X()[i] * conc;
            mix.jacobianRho(mix.T(), temp, sjac_exact);

            // Now compute the finite-differenced jacobian
            double eps = std::sqrt(Numerics::NumConst<double>::eps);
            double delta, prev;
            mix.netProductionRates(mix.T(), temp, species_values);
            for (int j = 0; j < mix.nSpecies(); ++j) {
                delta = std::max(temp[j] * eps, 1.0e-100);
                prev = temp[j];
                temp[j] += delta;
                mix.netProductionRates(mix.T(), temp, temp2);
                temp[j] = prev;
                for (int i = 0; i < mix.nSpecies(); ++i)
                    sjac_fd[i*mix.nSpecies()+j] = (temp2[i]-species_values[i])/
                        (delta*mix.speciesMw(j));
            }

            double sum1 = 0.0;
            double sum2 = 0.0;
            double x;
            for (int i = 0; i < mix.nSpecies(); ++i) {
                for (int j = 0; j < mix.nSpecies(); ++j) {
                    jac_file << setw(15)
                             //<< (sjac_fd[i*mix.nSpecies()+j] -
                             //       sjac_exact[i*mix.nSpecies()+j]) /
                             //       sjac_exact[i*mix.nSpecies()+j];
                             << sjac_exact[i*mix.nSpecies()+j];

                    // Compute the Frobenius norm of the difference matrix
                    //x = sjac_exact[i*mix.nSpecies()+j] - sjac_fd[i*mix.nSpecies()+j];
                    //sum1 += x*x;
                    //x = sjac_exact[i*mix.nSpecies()+j];
                    //sum2 += x*x;
                }
                jac_file << endl;
            }
            //jac_file << setw(10) << mix.T() << setw(15) << std::sqrt(sum1/sum2);
            jac_file << endl;*/
        }
    }

    // Clean up storage
    delete [] species_values;
    delete [] temp;
    delete [] sjac_exact;
    delete [] sjac_fd;
    delete [] temp2;
    delete [] reaction_values;

    return 0;
}
