
#include "cwrapper.h"
#include "mutation++.h"

#include <iostream>

using namespace std;
using namespace Mutation::Thermodynamics;

Mutation::Mixture* p_mix;
double* p_work_species;
double* p_work_element;
//RealVector work_species;

//==============================================================================
std::string char_to_string(char *str, int length)
{
    std::string string(str, length);
    string.erase(0, string.find_first_not_of(' '));
    string.erase(string.find_last_not_of(' ')+1);
    return string;
}

//==============================================================================
void string_to_char(std::string string, char* str, int length)
{
    for (int i = 0; i < length; ++i)
        if (i < string.length())
            str[i] = string[i];
        else
            str[i] = ' ';    
}

//==============================================================================
void NAME_MANGLE(initialize)(char* mixture, int mixture_length)
{
    //feenableexcept(FE_INVALID);
    p_mix = new Mutation::Mixture(char_to_string(mixture, mixture_length));
    p_work_species = new double [p_mix->nSpecies()];
    p_work_element = new double [p_mix->nElements()];
    //work_species = RealVector(p_mix->nSpecies());
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
int NAME_MANGLE(element_index)(char* element, int element_length)
{
    return p_mix->elementIndex(char_to_string(element, element_length)) + 1;
}

//==============================================================================
int NAME_MANGLE(species_index)(char* species, int species_length)
{
    return p_mix->speciesIndex(char_to_string(species, species_length)) + 1;
}

//==============================================================================
void NAME_MANGLE(species_name)(int* index, char* species, int species_length)
{
    string_to_char(p_mix->speciesName(*index-1), species, species_length);
}

//==============================================================================
double NAME_MANGLE(mixture_mw)()
{
    return p_mix->mixtureMw();
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
double NAME_MANGLE(number_density)()
{
    return p_mix->numberDensity();
}

//==============================================================================
double NAME_MANGLE(pressure)()
{
    return p_mix->P();
}

//==============================================================================
double NAME_MANGLE(density)()
{
    return p_mix->density();
}

//==============================================================================
void NAME_MANGLE(species_densities)(double* const rhoi)
{
    double rho = p_mix->density();
    for (int i = 0; i < p_mix->nSpecies(); ++i)
        rhoi[i] = rho * p_mix->Y()[i];
}

//==============================================================================
void NAME_MANGLE(equilibrate)(double* T, double* P)
{
    p_mix->equilibrate(*T, *P);
}

//==============================================================================
void NAME_MANGLE(set_state_t_rhoi)(double* T, double* rhoi)
{
    double rho = 0.0;
    for (int i = 0; i < p_mix->nSpecies(); ++i)
        rho += rhoi[i];
    
    double P = 0.0;
    for (int i = 0; i < p_mix->nSpecies(); ++i) {
        p_work_species[i] = rhoi[i] / rho;
        P += rhoi[i] / p_mix->speciesMw(i);
    }
    
    P *= Mutation::RU * (*T);
    p_mix->setStateTPY(T, &P, p_work_species);
}

//==============================================================================
void NAME_MANGLE(x)(double *const X)
{
    std::copy(p_mix->X(), p_mix->X()+p_mix->nSpecies(), X);
}

//==============================================================================
void NAME_MANGLE(y)(double *const Y)
{
    //std::copy(p_mix->Y(), p_mix->Y()+p_mix->nSpecies(), Y);
    for (int i = 0; i < p_mix->nSpecies(); ++i)
        Y[i] = p_mix->Y()[i];
}

//==============================================================================
void NAME_MANGLE(species_cp_mass)(double* const cp)
{
    p_mix->speciesCpOverR(cp);
    for (int i = 0; i < p_mix->nSpecies(); i++)
        cp[i] *= Mutation::RU / p_mix->speciesMw(i);
}

//==============================================================================
double NAME_MANGLE(mixture_frozen_cp_mass)()
{
    return p_mix->mixtureFrozenCpMass();
}

//==============================================================================
double NAME_MANGLE(mixture_frozen_cv_mass)()
{
    return p_mix->mixtureFrozenCvMass();
}

//==============================================================================
double NAME_MANGLE(mixture_frozen_gamma)()
{
    return p_mix->mixtureFrozenGamma();
}

//==============================================================================
void NAME_MANGLE(species_h_mass)(double *const h)
{
    p_mix->speciesHOverRT(h);
    double T = p_mix->T();
    for (int i = 0; i < p_mix->nSpecies(); i++)
        h[i] *= (Mutation::RU * T / p_mix->speciesMw(i));
}

//==============================================================================
double NAME_MANGLE(mixture_h_mass)()
{
    return p_mix->mixtureHMass();
}

//==============================================================================
double NAME_MANGLE(mixture_e_mass)()
{
    return p_mix->mixtureEnergyMass();
}

//==============================================================================
void NAME_MANGLE(net_production_rates)(double* const wdot)
{
    p_mix->netProductionRates(wdot);
}

//==============================================================================
int NAME_MANGLE(ncollision_pairs)()
{
    return p_mix->nCollisionPairs();
}


//==============================================================================
double NAME_MANGLE(viscosity)()
{
    return p_mix->viscosity();
}

//==============================================================================
double NAME_MANGLE(frozen_thermal_conductivity)()
{
    return p_mix->frozenThermalConductivity();
}

//==============================================================================
double NAME_MANGLE(equilibrium_thermal_conductivity)()
{
    return p_mix->equilibriumThermalConductivity();
}

//==============================================================================
double NAME_MANGLE(heavy_thermal_conductivity)()
{
    return p_mix->heavyThermalConductivity();
}

//==============================================================================
double NAME_MANGLE(electron_thermal_conductivity)()
{
    return p_mix->electronThermalConductivity();
}

//==============================================================================
double NAME_MANGLE(internal_thermal_conductivity)()
{
    return p_mix->internalThermalConductivity();
}

//==============================================================================
double NAME_MANGLE(reactive_thermal_conductivity)()
{
    return p_mix->reactiveThermalConductivity();
}

//==============================================================================
void NAME_MANGLE(average_diffusion_coeffs)(double *const p_Di)
{
    p_mix->averageDiffusionCoeffs(p_Di);
}

//==============================================================================
void NAME_MANGLE(stefan_maxwell)
    (const double* const p_dp, double* const p_V, double* const p_E)
{
    p_mix->stefanMaxwell(p_dp, p_V, *p_E);
}
    
//==============================================================================
double NAME_MANGLE(sigma)()
{
    return p_mix->sigma();
}

//==============================================================================
void NAME_MANGLE(surface_mass_balance)
    (const double *const p_Yke, const double *const p_Ykg, const double* const T,
     const double* const P, const double* const Bg, double* const Bc,
     double* const hw, double *const p_Xs)
{
    p_mix->surfaceMassBalance(p_Yke, p_Ykg, *T, *P, *Bg, *Bc, *hw, p_Xs);
}




