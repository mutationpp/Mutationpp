/**
 * @file cwrapper.cpp
 *
 * @brief Implements the C-wrapper to the Mutation++ library.
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
std::string char_to_string(F_STRING str, F_STRLEN length)
{
    std::string string(str, length);
    string.erase(0, string.find_first_not_of(' '));
    string.erase(string.find_last_not_of(' ')+1);
    return string;
}

//==============================================================================
void string_to_char(std::string string, F_STRING str, F_STRLEN length)
{
    for (int i = 0; i < length; ++i)
        if (i < string.length())
            str[i] = string[i];
        else
            str[i] = ' ';    
}

//==============================================================================
void NAME_MANGLE(initialize)(
    F_STRING mixture, F_STRING state_model, F_STRLEN mixture_length,
    F_STRLEN state_length)
{
    //feenableexcept(FE_INVALID);
    Mutation::MixtureOptions opts(char_to_string(mixture, mixture_length));
    opts.setStateModel(char_to_string(state_model, state_length));
    p_mix = new Mutation::Mixture(opts);
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
int NAME_MANGLE(n_mass_eqns)()
{
    return p_mix->nMassEqns();
}

//==============================================================================
int NAME_MANGLE(n_energy_eqns)()
{
    return p_mix->nEnergyEqns();
}

//==============================================================================
int NAME_MANGLE(element_index)(F_STRING element, F_STRLEN element_length)
{
    return p_mix->elementIndex(char_to_string(element, element_length)) + 1;
}

//==============================================================================
int NAME_MANGLE(species_index)(F_STRING species, F_STRLEN species_length)
{
    return p_mix->speciesIndex(char_to_string(species, species_length)) + 1;
}

//==============================================================================
void NAME_MANGLE(species_name)(int* index, F_STRING species, F_STRLEN species_length)
{
    string_to_char(p_mix->speciesName(*index-1), species, species_length);
}

//==============================================================================
double NAME_MANGLE(mixture_mw)()
{
    return p_mix->mixtureMw();
}

//==============================================================================
double NAME_MANGLE(mixture_t)()
{
    return p_mix->T();
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
void NAME_MANGLE(convert_xe_to_ye)(
    const double* const element_x, double* const element_y)
{
    p_mix->convert<XE_TO_YE>(element_x, element_y);
}
//==============================================================================
void NAME_MANGLE(convert_ye_to_xe)(
    const double* const element_y, double* const element_x)
{
    p_mix->convert<YE_TO_XE>(element_y, element_x);
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
void NAME_MANGLE(density_tpx)(double* T, double* P, const double* const X, double* rho)
{
    *rho = p_mix->density(*T, *P, X);
}

//==============================================================================
void NAME_MANGLE(species_densities)(double* const rhoi)
{
    double rho = p_mix->density();
    for (int i = 0; i < p_mix->nSpecies(); ++i)
        rhoi[i] = rho * p_mix->Y()[i];
}

//==============================================================================
void NAME_MANGLE(equilibrium_composition)(double* T, double* P, double* X)
{
    p_mix->equilibriumComposition(*T, *P, X);
}

//==============================================================================
void NAME_MANGLE(pyro_equilibrium_composition)(double* T, double* P, double* el, double* X)
{
    p_mix->equilibriumComposition(*T, *P, el, X);
}

//==============================================================================
void NAME_MANGLE(set_state)(double* v1, double* v2, int* vars)
{
    p_mix->setState(v1, v2, *vars);
}

//==============================================================================
void NAME_MANGLE(x)(double *const X)
{
    std::copy(p_mix->X(), p_mix->X()+p_mix->nSpecies(), X);
}

//==============================================================================
void NAME_MANGLE(y)(double *const Y)
{
    std::copy(p_mix->Y(), p_mix->Y()+p_mix->nSpecies(), Y);
}

//==============================================================================
void NAME_MANGLE(get_temperatures)(double* const T)
{
    p_mix->getTemperatures(T);
}

//==============================================================================
void NAME_MANGLE(species_cp_mass)(double* const cp)
{
    p_mix->getCpsMass(cp); 
}

//==============================================================================
void NAME_MANGLE(species_cv_mass)(double* const cv)
{
    p_mix->getCvsMass(cv); 
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
double NAME_MANGLE(mixture_frozen_sound_speed)()
{
    return p_mix->frozenSoundSpeed();
}

//==============================================================================
void NAME_MANGLE(species_e_mass)(double* const e)
{
    p_mix->getEnergiesMass(e);
}

//==============================================================================
void NAME_MANGLE(species_h_mass)(double *const h)
{
    p_mix->getEnthalpiesMass(h);
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
void NAME_MANGLE(frozen_thermal_conductivity)(double* const lambda)
{
    p_mix->frozenThermalConductivityVector(lambda); 
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
    (const double *const p_Ykc, const double *const p_Yke, const double *const p_Ykg, const double* const T,
     const double* const P, const double* const Bg, double* const Bc,
     double* const hw, double *const p_Xs)
{
    p_mix->surfaceMassBalance(p_Ykc, p_Yke, p_Ykg, *T, *P, *Bg, *Bc, *hw, p_Xs);
}

//==============================================================================
void NAME_MANGLE(source_energy_transfer)
     (double *const p_source_transfer)
{
     p_mix->energyTransferSource(p_source_transfer);
}
//==============================================================================
void NAME_MANGLE(surface_mass_balance_stag_line)
    (F_STRING mixture_name, double*  T, double*  P, double * p_Xg, double& ykf,F_STRLEN mixture_name_length)
{
   // string mixture_name = "Sio2";
   // string mixture_name = "Iron_meteor5";
    string name_name =char_to_string(mixture_name, mixture_name_length);
    string wall_name = "_wall";
    string total_mixture_name = name_name+wall_name;


    Mutation::Mixture mix_stagline (total_mixture_name);


    const int ne = mix_stagline.nElements();
    const int ns = mix_stagline.nSpecies();
    const int ng = mix_stagline.nGas();
    const int nc = mix_stagline.nCondensed();

    double* p_Ykc = new double [ne];
    double* p_Yke = new double [ne];
    double* p_Ykg = new double [ne];
    double* p_Xs = new double [ns];
    double* p_Ycond = new double [nc];


    mix_stagline.getComposition("BLedge", p_Yke, Composition::MASS);
    mix_stagline.getComposition("Gas", p_Ykg, Composition::MASS);
    mix_stagline.getComposition("Char" , p_Ykc, Composition::MASS);


    mix_stagline.surfaceMassBalance_CFD(p_Ykc, p_Yke, p_Ykg, *T, *P, p_Xs);

    for(int i = 0; i<ng;++i){
        p_Xg[i] = p_Xs[i];
    }

   /* for(int i = 0; i<ng;++i){
        cout<<"p_Xg("<<mix_stagline.speciesName(i)<<") :"<<p_Xg[i]<<endl;
    }*/

     mix_stagline.convert<X_TO_Y> (p_Xs,p_Xs);

     int spe =0;
     int ele =0;

     if (name_name == "Sio2"){
      spe = mix_stagline.speciesIndex("SiO2(L)");
      ele = mix_stagline.elementIndex("Si"); }

     else if(name_name == "Iron_meteor"){
          spe = mix_stagline.speciesIndex("Fe(L)");
          ele = mix_stagline.elementIndex("Fe");}


    ykf = (mix_stagline.elementMatrix()(spe,ele)*mix_stagline.atomicMass(ele)*p_Xs[spe])/ mix_stagline.speciesMw(spe);


    //ykf = 0.;p_Xs[spe

    delete [] p_Ycond;
    delete [] p_Yke;
    delete [] p_Ykg;
    delete [] p_Ykc;
    delete [] p_Xs;
}
//==============================================================================
void NAME_MANGLE(convert_ys_to_ye)(
    const double* species_y, double* elements_y)
{

    p_mix->convert<Y_TO_YE>(species_y, elements_y);
}

//==============================================================================
void NAME_MANGLE(number_of_elements)(int& element)
{
    element = p_mix->nElements();
}

//==============================================================================
     
     


