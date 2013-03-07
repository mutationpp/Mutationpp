/**
 * @file cwrapper.h
 *
 * This is the header file for the C wrapper to the Mutation++ library.  It
 * provides a subset of the functionality in C so that it may be used in Fortran
 * codes as well.
 */
 
#ifndef CWRAPPER_H
#define CWRAPPER_H

// Converts the base function name in mpp_function(_) where the (_) represents
// the name mangling that is performed in order to use the function in Fortran
#define NAME_MANGLE(__name__) mpp_##__name__##_

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
 * Below is an example Fortran snippit which illustrates how these functions may
 * be used.
 *
 * @code
 * program test
 *    use mutationpp ! Needed to use the wrapper
 *
 *    character(len=10) :: mixture
 *    integer :: ne, ns, nr
 *    real :: T, P, rho, n, cp, cv, h, mw, e
 *    real, dimension(:), allocatable :: element_x
 *    real, dimension(:), allocatable :: species_x
 *    real, dimension(:), allocatable :: species_y
 *
 *    ! Define the mixture to be used
 *    mixture = "air11"
 * 
 *    ! Initialize the Mutation++ library
 *    call mpp_initialize(mixture)
 *    ne = mpp_nelements()
 *    ns = mpp_nspecies()
 *    nr = mpp_nreactions()
 *    
 *    ! Allocate storage
 *    allocate(element_x(ne))
 *    allocate(species_x(ns))
 *    allocate(species_y(ns))
 *   
 *    ! Initialize elemental fractions for equilibrium calculation
 *    element = "N";  element_x(mpp_element_index(element)) = 0.79;
 *    element = "O";  element_x(mpp_element_index(element)) = 0.21;
 *    element = "e-"; element_x(mpp_element_index(element)) = 0.0;
 *
 *    ! Compute some properties at 1000K and 1atm
 *    T = 1000.0
 *    P = 101325.0
 *    call mpp_equilibrate_mole(T, P, element_x, species_x)
 *    call mpp_convert_x_to_y(species_x, species_y)        
 *    call mpp_number_density(T, P, n)
 *    call mpp_density(T, P, species_x, rho)
 *    call mpp_mixture_mw_mole(species_x, mw)
 *    call mpp_mixture_frozen_cp_mass(T, species_y, cp)
 *    call mpp_mixture_frozen_cv_mass(T, species_y, cv)
 *    call mpp_mixture_h_mass(T, species_y, h)
 *    call mpp_mixture_e_mass(T, rho, species_y, e)
 *
 *    ! Clean up the library memory (should be done once at the end)
 *    mpp_destroy()
 * end program test
 * @endcode
 * 
 * @{
 */

/**
 * Initialize the mutation++ library using a given mixture name.  This should be
 * called once before calling any other function.
 */
void NAME_MANGLE(initialize)(
    char* mixture, int mixture_length);

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
 * Returns the index of the element with the given name.
 */
int NAME_MANGLE(element_index)(
    char* element, int element_length);

/**
 * Returns the index of the species with the given name or -1 if the species
 * does not exist in the mixture.
 *
 * @param species - name of the species
 */
int NAME_MANGLE(species_index)(
    char* species, int species_length);

/**
 * Returns the name of the species with the given index.
 *
 * @param index   - index of the species (should be between 1 and nspecies
 *                  inclusive)
 * @param species - the name of the species on return
 */
void NAME_MANGLE(species_name)(
    int* index, char* species, int species_length);

/**
 * Returns the mixture molecular weight in kg/mol.
 */
double NAME_MANGLE(mixture_mw)();

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
 * Returns the density of the mixture given the mixture temperature and
 * pressure and species mole fractions.
 */
double NAME_MANGLE(density)();

/**
 * Returns the equilibrium composition of the mixture in species mass fractions
 * given the temperature, pressure, and elemental mass fractions.
 *
 * @param T         - mixture temperature
 * @param P         - mixture pressure
 * @param element_y - element densities or mass fractions
 * @param species_y - species mass fractions on return
 */
void NAME_MANGLE(equilibrate)(
    double* T, double* P);

/**
 * Returns the species specific heats at constant pressure in J/kg-K given the
 * mixture temperature.
 */
void NAME_MANGLE(species_cp_mass)(double* const cp);

/**
 * Returns the mixture specific heat at constant pressure in J/kg-K given the
 * mixture temperature and species mass fractions.
 */
double NAME_MANGLE(mixture_frozen_cp_mass)();

/**
 * Returns the mixture specific heat at constant volume in J/kg-K given the
 * mixture temperature and species mass fractions.
 */
double NAME_MANGLE(mixture_frozen_cv_mass)();

/**
 * Returns the mixture ratio of specific heats \f$C_p/C_v\f$ which is a unitless 
 * quantity.
 */
double NAME_MANGLE(mixture_frozen_gamma)();

/**
 * Returns the species enthalpies in J/kg.
 *
 * @param h - species enthalpies on return
 */
void NAME_MANGLE(species_h_mass)(double *const h);

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
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif // CWRAPPER_H
