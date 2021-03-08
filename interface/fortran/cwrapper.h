/**
 * @file cwrapper.h
 *
 * @brief This is the header file for the C-wrapper to the Mutation++ library.
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

#ifndef FORTRAN_CWRAPPER_H
#define FORTRAN_CWRAPPER_H

// Converts the base function name in mpp_function(_) where the (_) represents
// the name mangling that is performed in order to use the function in Fortran
#define NAME_MANGLE(__name__) mpp_##__name__##_

#define F_STRING char*
#define F_STRLEN long int

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup FortranWrapper Fortran Wrapper
 *
 * This wrapper provides a set of C functions which wrap the Mutation++ library
 * and can be called inside of a Fortran code by simply linking the 
 * libmutation_fortran.a library with the fortran code.  There are two things to
 * note about the function declarations shown below.
 *
 * - The functions are named with an underscored appended at the end of each 
 * function. This is done in order to match the name mangling performed by the 
 * Fortran compiler.  If your Fortran compiler uses a different name mangling
 * procedure, then you will have to change this in cwrapper.h by altering the
 * NAME_MANGLE macro at the top of the file.
 *
 * - The function declarations include length parameters for any string that is
 * used as input or output to the functions.  These do not need to be included
 * when calling the function in a Fortran code.
 *
 * See \ref wrapper_test.f90 for an example of a code which uses this wrapper.
 *
 * @{
 */

/**
 * Initialize the mutation++ library using a given mixture name.  This should be
 * called once before calling any other function.
 */
void NAME_MANGLE(initialize)(
    F_STRING mixture, F_STRING state_model, F_STRLEN mixture_length,
    F_STRLEN state_length);

/**
 * Deallocates all data associated with the mutation++ library.  Should be
 * called once after all other functions have been called.
 */
void NAME_MANGLE(destroy)();

/**
 * Returns the number of elements in the mixture.
 */
int NAME_MANGLE(nelements)();

/**
 * Returns the number of species in the mixture.
 */
int NAME_MANGLE(nspecies)();

/**
 * Returns the number of reactions in the mixture.
 */
int NAME_MANGLE(nreactions)();

/**
 * Returns the number of mass equations associated with the state model.
 */
int NAME_MANGLE(n_mass_eqns)();

/**
 * Returns the number of energy equations associated with the state model.
 */
int NAME_MANGLE(n_energy_eqns)();

/**
 * Returns the index of the element with the given name.
 */
int NAME_MANGLE(element_index)(
    F_STRING element, F_STRLEN element_length);

/**
 * Returns the index of the species with the given name or -1 if the species
 * does not exist in the mixture.
 *
 * @param species - name of the species
 */
int NAME_MANGLE(species_index)(
    F_STRING species, F_STRLEN species_length);

/**
 * Returns the name of the species with the given index.
 *
 * @param index   - index of the species (should be between 1 and nspecies
 *                  inclusive)
 * @param species - the name of the species on return
 */
void NAME_MANGLE(species_name)(
    int* index, F_STRING species, F_STRLEN species_length);

/**
 * Returns the mixture molecular weight in kg/mol.
 */
double NAME_MANGLE(mixture_mw)();

/**
 * Returns the mixture translational temperature in K.
 */
double NAME_MANGLE(mixture_t)();

/**
 * Returns the array of species molecular weights in kg/mol.
 *
 * @param mw - the species molecular weights on return
 */
void NAME_MANGLE(species_mw)(double* const mw);

/**
 * Converts the species mole fractions to species mass fractions.
 *
 * @param species_x - species mole fractions
 * @param species_y - species mass fractions on return
 */
void NAME_MANGLE(convert_x_to_y)(
    const double* species_x, double* species_y);

/**
 * Converts the element mole fractions to element mass fractions.
 *
 * @param element_x - element mole fractions
 * @param element_y - element mass fractions on return
 */
void NAME_MANGLE(convert_xe_to_ye)(
    const double* element_x, double* element_y);

/**
 * Converts the species densities to species mole fractions.
 *
 * @param species_rho - species mass densities
 * @param species_x   - species mole fractions on return
 */
void NAME_MANGLE(convert_rho_to_x)(
    const double* species_rho, double* species_x);

/**
 * Computes the number of moles of each element that exist in a mixture with
 * the given set of species moles.
 *
 * @param species_N - species moles or mole fractions
 * @param element_N - element moles or mole fractions
 */
void NAME_MANGLE(element_moles)(
    const double *const species_N, double *const element_N);

/**
 * Returns the mole fractions for the current mixture state.
 */
void NAME_MANGLE(x)(double* const X);

/**
 * Returns the mass fractions for the current mixture state.
 */
void NAME_MANGLE(y)(double* const Y);

/**
 * Fills the temperature vector with the current mixture state.
 */
void NAME_MANGLE(get_temperatures)(double* const T);

/**
 * Returns the number density of the mixture given the mixture temperature
 * and pressure.
 */
double NAME_MANGLE(number_density)();

/**
 * Returns the pressure of the mixture given the mixture temperature and
 * density and species mass fractions.
 */
double NAME_MANGLE(pressure)();

/**
 * Returns the density of the mixture.
 */
double NAME_MANGLE(density)();

/**
 * Returns the density of the mixture given the mixture temperature and
 * pressure and species mole fractions.
 */
void NAME_MANGLE(density_tpx)(double* T, double* P, const double* const X, double* rho);

/**
 * Returns the current species densities.
 */
void NAME_MANGLE(species_densities)(
    double* const rhoi);

/**
 * Returns the equilibrium composition of the mixture in species mole fractions
 * given the temperature and pressure.
 *
 * @param T - mixture temperature
 * @param P - mixture pressure
 * @param X - mole fractions
 */
void NAME_MANGLE(equilibrium_composition)(double* T, double* P, double* X);
void NAME_MANGLE(pyro_equilibrium_composition)(double* T, double* P, double* el, double* X);

/**
 * Sets the current state of the mixture with the appropriate variable set.
 */
void NAME_MANGLE(set_state)(double* v1, double* v2, int* vars);

/**
 * Returns the species specific heats at constant pressure in J/kg-K per species
 * per each energy mode.
 */
void NAME_MANGLE(species_cp_mass)(double* const cp);

/**
 * Returns the species specific heats at constant volume in J/kg-K per species
 * per each energy mode.
 */
void NAME_MANGLE(species_cv_mass)(double* const cv);

/**
 * Returns the mixture specific heat at constant pressure in J/kg-K.
 */
double NAME_MANGLE(mixture_frozen_cp_mass)();

/**
 * Returns the mixture specific heat at constant volume in J/kg-K.
 */
double NAME_MANGLE(mixture_frozen_cv_mass)();

/**
 * Returns the mixture ratio of specific heats \f$C_p/C_v\f$ which is a unitless 
 * quantity.
 */
double NAME_MANGLE(mixture_frozen_gamma)();

/**
 * Returns the mixture frozen sound speed in m/s.
 */
double NAME_MANGLE(mixture_frozen_sound_speed)();

/**
 * Returns the species energies (total + internal if multi temperature) in J/kg.
 *
 * @param e - species energies on return
 */
void NAME_MANGLE(species_e_mass)(double *const e);

/**
 * Returns the species enthalpies (total + internal if multi temperature) in J/kg.
 *
 * @param h - species enthalpies on return
 */
void NAME_MANGLE(species_h_mass)(double *const h);

/**
 * Returns the species entropies (total + internal if multi temperature) in J/kg-K.
 *
 * @param s - species entropies on return
 */
void NAME_MANGLE(species_s_mass)(double *const s);

/**
 * Returns the mixture enthalpy in J/kg given the mixture temperature and
 * species mass fractions.
 */
double NAME_MANGLE(mixture_h_mass)();

/**
 * Return the mixture total energy in J/kg given temperature, density, and mass
 * fractions.
 */
double NAME_MANGLE(mixture_e_mass)();

/**
 * Fills the vector wdot with the net species production rates due to the
 * chemical reactions.
 * \f[
 * \dot{\omega}_i = M_{w,i} \sum_j \left( \nu_{ij}^" - \nu_{ij}^{'} \right) 
 *    \left[ k_{f,j} \prod_i C_i^{\nu_{ij}^{'}} - k_{b,j} \prod_i 
 *    C_i^{\nu_{ij}^"} \right] \Theta_{TB}
 * \f]
 *
 * @param wdot  on return, the species production rates in kg/m^3-s
 */
void NAME_MANGLE(net_production_rates)(double* const wdot);

/**
 * Fills the matrix j with the species production rate jacobian matrix
 * \f[
 * J_{ij} = \frac{\partial \dot{\omega}_i}{\partial \rho_j}
 * \f]
 * The Jacobian matrix should be sized to be at least the square of the
 * number of species, ns.  Access the jacobian using row-major ordering (ie:
 * J_{ij} = p_jac[i*ns + j]).
 *
 * @param p_jac  - on return, the jacobian matrix \f$J_{ij}\f$
 */
void NAME_MANGLE(species_jacobian_rho)(double* const j);

/**
 * Returns the total number of collision pairs accounted for in the 
 * collision database.
 */
int NAME_MANGLE(ncollision_pairs)();

/**
 * Returns the mixture viscosity.
 */
double NAME_MANGLE(viscosity)();

/**
 * Returns the mixture thermal conductivity for a frozen mixture.
 */
void NAME_MANGLE(frozen_thermal_conductivity)(double* const lambda);

/**
 * Returns the heavy thermal diffusion ratios for each species.
 */
void NAME_MANGLE(heavy_thermal_diffusion_ratios)(double* const pk);

/**
 * Returns the mixture thermal conductivity for a mixture in thermochemical
 * equilibrium.
 */
double NAME_MANGLE(equilibrium_thermal_conductivity)();

/**
 * Returns the heavy particle translational thermal conductivity using the 
 * set algorithm.
 */
double NAME_MANGLE(heavy_thermal_conductivity)();

/**
 * Returns the electron translational thermal conductivity.
 */
double NAME_MANGLE(electron_thermal_conductivity)();

/**
 * Returns the internal energy thermal conductivity using Euken's formulas.
 */
double NAME_MANGLE(internal_thermal_conductivity)(double T);

/**
 * Returns the reactive thermal conductivity which accounts for reactions
 * for mixtures in thermochemical equilibrium.
 */
double NAME_MANGLE(reactive_thermal_conductivity)();


/**
 * Returns mixture averaged species diffusion coefficients which are defined
 * as \f[ D_i = \frac{1 - X_i}{\sum_j X_j / \mathcal{D}_{ij}} \f] where
 * \f$\mathcal{D}_{ij}\f$ are the binary diffusion coefficients for each
 * species pair.
 */
void NAME_MANGLE(average_diffusion_coeffs)(double *const p_Di);

/**
 * Computes the species diffusion velocities and ambipolar electric field
 * using the Ramshaw approximation of the generalized Stefan-Maxwell
 * equations and the supplied modified driving forces
 *
 * \f[ d^{'}_i = \frac{\nabla p_i}{n k_B T_h} - \frac{y_i p}{n k_B T_h}
 *     \nabla \ln p + k^h_{Ti} \nabla \ln T_h + k^e_{Ti} \frac{T_h}{T_e}
 *     \nabla \ln T_e \f]
 *
 * @param p_dp - the vector of modified driving forces \f$ d^{'}_i \f$
 * @param p_V  - on return, the vector of diffusion velocities
 * @param E    - on return, the ambipolar electric field
 */
void NAME_MANGLE(stefan_maxwell)
    (const double* const p_dp, double* const p_V, double* const p_E);

void NAME_MANGLE(diffusion_matrix)(double* const p_Dij);

/**
 * Returns the electric conductivity.
 */
double NAME_MANGLE(sigma)();

/**
 * Solves the surface mass balance at an ablating surface.
 */
void NAME_MANGLE(surface_mass_balance)
    (const double *const p_Yke, const double *const p_Ykg, const double* const T,
     const double* const P, const double* const Bg, double* const Bc,
     double* const hw, double *const p_Xs);

/**
 * Returns the pointer to the energy transfer between the internal
 * energy modes.
 */

void NAME_MANGLE(source_energy_transfer)
    (double* const p_source_transfer);

/**
 * Sets the surface state using appropriate variable set.
 */
void NAME_MANGLE(set_surface_state)(double* v1, double* v2, int* vars);

/**
 * Returns the chemical production rates due to gas-surface interaction
 * reactions.
 */
void NAME_MANGLE(surface_production_rates)(double* v1);

/**
 * Sets the mole fractions at a distance dx from the surface needed for
 * the surface mass balance.
 */
void NAME_MANGLE(set_diffusion_model)( double* xi_edge, double* dx );

/**
 * Sets the temperature at a distance dx from the surface needed for
 * the surface energy balance.
 */
void NAME_MANGLE(set_cond_heat_flux)(double* T_edge, double* dx);

/**
 * Requests the solution of the surface balances at a surface using the
 * gas-surface interaction module of mutation++.
 */
void NAME_MANGLE(solve_surface_balance)();

/**
 * Returns the current surface state with respect to partial densities and
 * surface temperatures.
 */
void NAME_MANGLE(get_surface_state)(double* v1, double* v2, int* vars);

/**
 * Returns the mass blowing rate due to gas-surface interaction phenomena.
 */
void NAME_MANGLE(mass_blowing_rate)(double* mdot);

/**
 * Converts the element mass fractions to element mole fractions.
 *
 * @param element_y - element mass fractions
 * @param element_x - element mole fractions on return
 */
void NAME_MANGLE(convert_ye_to_xe)(
        const double* element_y, double* element_x);

/**
 * Converts the species mass fractions to element mass fractions.
 *
 * @param species_y - species mass fractions
 * @param elements_y - element mass fractions
 */

void NAME_MANGLE(convert_ys_to_ye)(
        const double* species_y, double* elements_y);

#ifdef __cplusplus
}
#endif

#endif // CWRAPPER_H
