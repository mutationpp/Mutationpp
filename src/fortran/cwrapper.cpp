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

#include <Eigen/Dense>

#include <iostream>

using namespace std;
using namespace Mutation::Thermodynamics;

Mutation::Mixture* p_mix;
double* p_work_species;
double* p_work_element;
double* p_energy_species_tot;
double* p_energy_species_trans;
double* p_energy_species_rot;
double* p_energy_species_vib;
double* p_energy_species_ele;
double* p_energy_species_form;

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
    p_energy_species_tot = new double [p_mix->nSpecies()];
    p_energy_species_trans = new double [p_mix->nSpecies()];
    p_energy_species_rot = new double [p_mix->nSpecies()];
    p_energy_species_vib = new double [p_mix->nSpecies()];
    p_energy_species_ele = new double [p_mix->nSpecies()];
    p_energy_species_form = new double [p_mix->nSpecies()];
    //work_species = RealVector(p_mix->nSpecies());
}

//==============================================================================
void NAME_MANGLE(destroy)()
{
    delete p_mix;
    delete [] p_work_species;
    delete [] p_work_element;
    
    delete [] p_energy_species_tot;
    delete [] p_energy_species_trans;
    delete [] p_energy_species_rot;
    delete [] p_energy_species_vib;
    delete [] p_energy_species_ele;
    delete [] p_energy_species_form;

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
void NAME_MANGLE(convert_y_to_x)(
    const double *const species_y, double *const species_x)
{
    p_mix->convert<Y_TO_X>(species_y, species_x);
}

//==============================================================================
void NAME_MANGLE(convert_xe_to_ye)(
    const double* const element_x, double* const element_y)
{
    p_mix->convert<XE_TO_YE>(element_x, element_y);
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
double NAME_MANGLE(internal_thermal_conductivity)(double T)
{
    return p_mix->internalThermalConductivity(T);
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

//==============================================================================
void NAME_MANGLE(source_energy_transfer)
     (double *const p_source_transfer)
{
     p_mix->energyTransferSource(p_source_transfer);
}

//==============================================================================
void NAME_MANGLE(set_wall_state)(double* v1, double* v2, int* vars)
{
     p_mix->setWallState(v1, v2, *vars);
}

//==============================================================================
void NAME_MANGLE(wall_production_rates)(double* const v1)
{
     p_mix->surfaceProductionRates(v1);
}

//==============================================================================
void NAME_MANGLE(set_diffusion_model)(double* rhoi_edge, double* dx)
{
     p_mix->setDiffusionModel(rhoi_edge, *dx);
}

//==============================================================================
void NAME_MANGLE(solve_surface_balance)()
{
     p_mix->solveSurfaceBalance();
}

//==============================================================================
void NAME_MANGLE(get_wall_state)(double* v1, double* v2, int* vars)
{
     p_mix->getWallState( v1, v2, *vars );
}
//==============================================================================
double NAME_MANGLE(r_univ)()
{
     return Mutation::RU;
}
//==============================================================================
double NAME_MANGLE(pi)()
{
     return Mutation::PI;
}
//==============================================================================
double NAME_MANGLE(kb)()
{
     return Mutation::KB;
}
//==============================================================================
void NAME_MANGLE(species_h_ref_form_mass)(double* h_tot_ref_mass, double* h_form_mass)
{
     const double T_ref = 298.15;
     p_mix->speciesHOverRT(T_ref, T_ref, T_ref, T_ref, T_ref, h_tot_ref_mass, NULL, NULL, NULL, NULL, h_form_mass);
     for (int i_sp = 0; i_sp<p_mix->nSpecies(); ++i_sp){
         h_tot_ref_mass[i_sp] *= ((Mutation::RU/p_mix->speciesMw(i_sp)) * T_ref);
         h_form_mass[i_sp] *= ((Mutation::RU/p_mix->speciesMw(i_sp)) * T_ref);
     }
}

//==============================================================================
void NAME_MANGLE(species_i_cv_mass_cosmic)(const int& i_fortran, double* T, double* cv)
{
    int i = i_fortran-1;

    p_mix->speciesCvOverR(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele);
    cv[0] = p_energy_species_tot[i] * (Mutation::RU / p_mix->speciesMw(i));
    cv[1] = p_energy_species_trans[i] * (Mutation::RU / p_mix->speciesMw(i));
    cv[2] = p_energy_species_rot[i] * (Mutation::RU / p_mix->speciesMw(i));
    cv[3] = p_energy_species_vib[i] * (Mutation::RU / p_mix->speciesMw(i));
    cv[4] = p_energy_species_ele[i] * (Mutation::RU / p_mix->speciesMw(i));
}

//==============================================================================
void NAME_MANGLE(species_i_cp_mass_cosmic)(const int& i_fortran, double* T, double* cp)
{
    int i = i_fortran-1;

    p_mix->speciesCpOverR(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele);
    cp[0] = p_energy_species_tot[i] * (Mutation::RU / p_mix->speciesMw(i));
    cp[1] = p_energy_species_trans[i] * (Mutation::RU / p_mix->speciesMw(i));
    cp[2] = p_energy_species_rot[i] * (Mutation::RU / p_mix->speciesMw(i));
    cp[3] = p_energy_species_vib[i] * (Mutation::RU / p_mix->speciesMw(i));
    cp[4] = p_energy_species_ele[i] * (Mutation::RU / p_mix->speciesMw(i));
}

//==============================================================================
void NAME_MANGLE(species_i_cv_mole_cosmic)(const int& i_fortran, double* T, double* cv)
{
    int i = i_fortran-1;

    p_mix->speciesCvOverR(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele);
    cv[0] = p_energy_species_tot[i] * (Mutation::RU);
    cv[1] = p_energy_species_trans[i] * (Mutation::RU);
    cv[2] = p_energy_species_rot[i] * (Mutation::RU);
    cv[3] = p_energy_species_vib[i] * (Mutation::RU);
    cv[4] = p_energy_species_ele[i] * (Mutation::RU);
}

//==============================================================================
void NAME_MANGLE(species_i_cp_mole_cosmic)(const int& i_fortran, double* T, double* cp)
{
    int i = i_fortran-1;

    p_mix->speciesCpOverR(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele);
    cp[0] = p_energy_species_tot[i] * (Mutation::RU);
    cp[1] = p_energy_species_trans[i] * (Mutation::RU);
    cp[2] = p_energy_species_rot[i] * (Mutation::RU);
    cp[3] = p_energy_species_vib[i] * (Mutation::RU);
    cp[4] = p_energy_species_ele[i] * (Mutation::RU);
}

//==============================================================================
void NAME_MANGLE(species_i_energy_mole_cosmic)(const int& i_fortran, double* T, double* energy)
{
    int i = i_fortran-1;

    p_mix->speciesHOverRT(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele, p_energy_species_form);
    energy[0] = (p_energy_species_tot[i]-1.0) * (Mutation::RU) * *T;
    energy[1] = (p_energy_species_trans[i]-1.0) * (Mutation::RU) * *T;
    energy[2] = p_energy_species_rot[i] * (Mutation::RU) * *T;
    energy[3] = p_energy_species_vib[i] * (Mutation::RU) * *T;
    energy[4] = p_energy_species_ele[i] * (Mutation::RU) * *T;
    energy[7] = p_energy_species_form[i] * (Mutation::RU) * *T;
}
 
//==============================================================================
void NAME_MANGLE(species_i_energy_mass_cosmic)(const int& i_fortran, double* T, double* energy)
{
    int i = i_fortran-1;

    p_mix->speciesHOverRT(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele, p_energy_species_form);
    energy[0] = (p_energy_species_tot[i]-1) * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    energy[1] = (p_energy_species_trans[i]-1.0) * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    energy[2] = p_energy_species_rot[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    energy[3] = p_energy_species_vib[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    energy[4] = p_energy_species_ele[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    energy[7] = p_energy_species_form[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
}

//==============================================================================
void NAME_MANGLE(species_i_enthalpy_mole_cosmic)(const int& i_fortran, double* T, double* enthalpy)
{
    int i = i_fortran-1;

    p_mix->speciesHOverRT(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele, p_energy_species_form);
    enthalpy[0] = (p_energy_species_tot[i]) * (Mutation::RU) * *T;
    enthalpy[1] = (p_energy_species_trans[i]) * (Mutation::RU) * *T;
    enthalpy[2] = p_energy_species_rot[i] * (Mutation::RU) * *T;
    enthalpy[3] = p_energy_species_vib[i] * (Mutation::RU) * *T;
    enthalpy[4] = p_energy_species_ele[i] * (Mutation::RU) * *T;
    enthalpy[7] = p_energy_species_form[i] * (Mutation::RU) * *T;
}
 
//==============================================================================
void NAME_MANGLE(species_i_enthalpy_mass_cosmic)(const int& i_fortran, double* T, double* enthalpy)
{
    int i = i_fortran-1;

    p_mix->speciesHOverRT(*T, *T, *T, *T, *T, 
                          p_energy_species_tot, p_energy_species_trans, p_energy_species_rot, 
                          p_energy_species_vib, p_energy_species_ele, p_energy_species_form);
    enthalpy[0] = (p_energy_species_tot[i]) * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    enthalpy[1] = (p_energy_species_trans[i]) * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    enthalpy[2] = p_energy_species_rot[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    enthalpy[3] = p_energy_species_vib[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    enthalpy[4] = p_energy_species_ele[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
    enthalpy[7] = p_energy_species_form[i] * (Mutation::RU / p_mix->speciesMw(i)) * *T;
}

//==============================================================================
void NAME_MANGLE(binary_diff_coeff)(double* Dbin)
{
    Eigen::MatrixXd m_Dbin = p_mix->binaryDiffusionCoefficients();

    for (int j = 0; j < p_mix->nSpecies(); j++){
        for (int i = 0; i < p_mix->nSpecies(); i++){
            if (i>j){
            Dbin[j*p_mix->nSpecies()+i] = m_Dbin(i,j);
            } else {
            Dbin[j*p_mix->nSpecies()+i] = m_Dbin(j,i);
            }
        }
    }
}
//==============================================================================
void NAME_MANGLE(jac_w_dot)(double* jac_w_dot)
{
     p_mix->jacobianRho(jac_w_dot);
}

