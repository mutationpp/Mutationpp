#include "mutation++.h"
#include "cwrapper.h"

#include <iostream>
using namespace std;

//#define _GNU_SOURCE
//#include <fenv.h>


using namespace Numerics;

Mixture* p_mix;
double* p_work_species;
double* p_work_element;
RealVector work_species;

//==============================================================================
std::string char_to_string(char *str, size_t length)
{
    std::string string(str, length);
    string.erase(0, string.find_first_not_of(' '));
    string.erase(string.find_last_not_of(' ')+1);
    return string;
}

//==============================================================================
void string_to_char(std::string string, char* str, size_t length)
{
    for (int i = 0; i < length; ++i)
        if (i < string.length())
            str[i] = string[i];
        else
            str[i] = ' ';    
}

//==============================================================================
void NAME_MANGLE(initialize)(char* mixture, size_t mixture_length)
{
    //feenableexcept(FE_INVALID);
    p_mix = new Mixture(char_to_string(mixture, mixture_length));
    p_work_species = new double [p_mix->nSpecies()];
    p_work_element = new double [p_mix->nElements()];
    work_species = RealVector(p_mix->nSpecies());
}

//==============================================================================
void NAME_MANGLE(destroy)()
{
    delete p_mix;
    delete [] p_work_species;
    delete [] p_work_element;
    
    p_mix = NULL;
    p_work_species = NULL;
    p_work_element = NULL;
}

//==============================================================================
int NAME_MANGLE(nelements)()
{
    return p_mix->nElements();
}

//==============================================================================
int NAME_MANGLE(nspecies)()
{
    return p_mix->nSpecies();
}

//==============================================================================
int NAME_MANGLE(nreactions)()
{
    return p_mix->nReactions();
}

//==============================================================================
int NAME_MANGLE(element_index)(char* element, size_t element_length)
{
    return p_mix->elementIndex(char_to_string(element, element_length)) + 1;
}

//==============================================================================
int NAME_MANGLE(species_index)(char* species, size_t species_length)
{
    return p_mix->speciesIndex(char_to_string(species, species_length)) + 1;
}

//==============================================================================
void NAME_MANGLE(species_name)(int* index, char* species, size_t species_length)
{
    string_to_char(p_mix->speciesName(*index-1), species, species_length);
}

//==============================================================================
void NAME_MANGLE(mixture_mw_mass)(
    const double *const Y, double *const Mw)
{
    *Mw = p_mix->mixtureMwMass(Y);
}

//==============================================================================
void NAME_MANGLE(mixture_mw_mole)(
    const double *const X, double *const Mw)
{
    *Mw = p_mix->mixtureMwMole(X);
}

//==============================================================================
void NAME_MANGLE(species_mw)(double* const mw)
{
    for (int i = 0; i < p_mix->nSpecies(); ++i)
        mw[i] = p_mix->speciesMw(i);
}

//==============================================================================
void NAME_MANGLE(convert_x_to_y)(
    const double *const species_x, double *const species_y)
{
    p_mix->convert<X_TO_Y>(species_x, species_y);
}

//==============================================================================
void NAME_MANGLE(convert_rho_to_x)(
    const double* species_rho, double* species_x)
{
    p_mix->convert<RHO_TO_X>(species_rho, species_x);
}

//==============================================================================
void NAME_MANGLE(element_moles)(
    const double *const species_N, double *const element_N)
{
    p_mix->elementMoles(species_N, element_N);
}

//==============================================================================
void NAME_MANGLE(number_density)(
    const double *const T, const double *const P, double *const n)
{
    *n = p_mix->numberDensity(*T, *P);
}

//==============================================================================
void NAME_MANGLE(pressure)(
    const double *const T, const double *const rho, const double *const Y, 
    double *const P)
{
    *P = p_mix->pressure(*T, *rho, Y);
}

//==============================================================================
void NAME_MANGLE(density)(
    const double *const T, const double *const P, const double *const X,
    double *const rho)
{
    *rho = p_mix->density(*T, *P, X);
}

//==============================================================================
void NAME_MANGLE(equilibrate_mass)(
    double* T, double* P, double* element_y, double* species_y)
{
    p_mix->convert<YE_TO_XE>(element_y, p_work_element);
    p_mix->equilibrate(*T, *P, p_work_element, species_y);
    p_mix->convert<X_TO_Y>(species_y, species_y);
}

//==============================================================================
void NAME_MANGLE(equilibrate_mole)(
    double* T, double* P, double* element_x, double* species_x)
{
    p_mix->equilibrate(*T, *P, element_x, species_x);
}

//==============================================================================
void NAME_MANGLE(species_cp_mass)(
    const double* const T, double* const cp)
{
    p_mix->speciesCpOverR(*T, cp);
    for (int i = 0; i < p_mix->nSpecies(); i++)
        cp[i] *= RU / p_mix->speciesMw(i);
}

//==============================================================================
void NAME_MANGLE(mixture_frozen_cp_mass)(
    const double *const T, const double *const y, double *const cp)
{
    *cp = p_mix->mixtureFrozenCpMass(*T, y);
}

//==============================================================================
void NAME_MANGLE(mixture_frozen_cv_mass)(
    const double *const T, const double *const y, double *const cv)
{
    *cv = p_mix->mixtureFrozenCvMass(*T, y);
}

//==============================================================================
void NAME_MANGLE(mixture_frozen_gamma)(
    const double *const T, const double *const Y, double *const gamma)
{
    *gamma = p_mix->mixtureFrozenGamma(*T, Y);
}

//==============================================================================
void NAME_MANGLE(species_h_mass)(
    const double *const T, double *const h)
{
    p_mix->speciesHOverRT(*T, h);
    for (int i = 0; i < p_mix->nSpecies(); i++)
        h[i] *= (RU * *T / p_mix->speciesMw(i));
}

//==============================================================================
void NAME_MANGLE(mixture_h_mass)(
    const double *const T, const double *const y, double *const h)
{
    *h = p_mix->mixtureHMass(*T, y);
}

//==============================================================================
void NAME_MANGLE(mixture_e_mass)(
    const double *const T, const double *const rho, const double *const Y, 
    double *const e)
{
    *e = p_mix->mixtureHMass(*T, Y) - p_mix->pressure(*T, *rho, Y) / *rho;
}

//==============================================================================
void NAME_MANGLE(net_production_rates)(
    const double* const T, const double* const conc, double* const wdot)
{
    p_mix->netProductionRates(*T, conc, wdot);
}

