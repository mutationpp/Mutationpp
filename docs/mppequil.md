<a id="top"></a>

# Mutation++ Equilibrium (mppequil)
This tool essentially wraps the Mutation++ library to provide equilibrium properties for mixtures over a range of temperatures and/or pressures.  This is similar to tools like [NASA's CEA](https://www.grc.nasa.gov/www/CEAWeb/).  

## Usage
```bash
mppequil [OPTIONS] mixture
Compute equilibrium properties for mixture over a set of temperatures and pressures using the Mutation++ library.

    -h, --help          prints this help message
        --no-header     no table header will be printed
    -T                  temperature range in K "T1:dT:T2" or simply T (default = 300:100:15000 K)
    -P                  pressure range in Pa "P1:dP:P2" or simply P (default = 1 atm)
    -B                  magnitude of the magnetic field in teslas (default = 0 T)
    -m                  list of mixture values to output (see below)
    -s                  list of species values to output (see below)
    -r                  list of reaction values to output (see below)
    -o                  list of other values to output (see below)
        --species-list  instead of mixture name, use this to list species in mixture
        --elem-x        set elemental mole fractions (ex: N:0.8,O:0.2)
        --elem-comp     set elemental composition with a name from the mixture file
        --thermo-db     overrides thermodynamic database type (NASA-7, NASA-9, RRHO)
        --scientific    outputs in scientific format with given precision

Mixture values (example format: "1-3,7,9-11"):
    0 : Th        [K]         heavy particle temperature
    1 : P         [Pa]        pressure
    2 : B         [T]         magnitude of magnetic field
    3 : rho       [kg/m^3]    density
    4 : nd        [1/m^3]     number density
    5 : Mw        [kg/mol]    molecular weight
    6 : Cp_eq     [J/mol-K]   equilibrium specific heat at constant pressure
    7 : H         [J/mol]     mixture enthalpy
    8 : S         [J/mol-K]   entropy
    9 : Cp_eq     [J/kg-K]    equilibrium specific heat at constant pressure
    10: H         [J/kg]      mixture enthalpy
    11: H-H0      [J/kg]      mixture enthalpy minus the enthalpy at 0K
    12: S         [J/kg-K]    entropy
    13: Cv_eq     [J/kg-K]    equilibrium specific heat at constant volume
    14: Cp        [J/mol-K]   frozen specific heat at constant pressure
    15: Cv        [J/mol-K]   frozen specific heat at constant volume
    16: Cp        [J/kg-K]    frozen specific heat at constant pressure
    17: Cv        [J/kg-K]    frozen specific heat at constant volume
    18: gam_eq    [-]         equilibrium ratio of specific heats
    19: gamma     [-]         frozen ratio of specific heat
    20: Ht        [J/mol]     translational enthalpy
    21: Hr        [J/mol]     rotational enthalpy
    22: Hv        [J/mol]     vibrational enthalpy
    23: Hel       [J/mol]     electronic enthalpy
    24: Hf        [J/mol]     formation enthalpy
    25: Ht        [J/kg]      translational enthalpy
    26: Hr        [J/kg]      rotational enthalpy
    27: Hv        [J/kg]      vibrational enthalpy
    28: Hel       [J/kg]      electronic enthalpy
    29: Hf        [J/kg]      formation enthalpy
    30: e         [J/mol]     mixture energy
    31: e         [J/kg]      mixture energy
    32: mu        [Pa-s]      dynamic viscosity
    33: lambda    [W/m-K]     mixture equilibrium thermal conductivity
    34: lam_reac  [W/m-K]     reactive thermal conductivity
    35: lam_bb    [W/m-K]     Butler-Brokaw reactive thermal conductivity
    36: lam_soret [W/m-K]     Soret thermal conductivity
    37: lam_int   [W/m-K]     internal energy thermal conductivity
    38: lam_h     [W/m-K]     heavy particle translational thermal conductivity
    39: lam_e     [W/m-K]     electron translational thermal conductivity
    40: sigma     [S/m]       electric conductivity (B=0)
    41: a_f       [m/s]       frozen speed of sound
    42: a_eq      [m/s]       equilibrium speed of sound
    43: Eam       [V/K]       ambipolar electric field (SM Ramshaw)
    44: drho/dP   [kg/J]      equilibrium density derivative w.r.t pressure

Species values (same format as mixture values):
    0 : X         [-]         mole fractions
    1 : dX/dT     [1/K]       partial of mole fraction w.r.t. temperature
    2 : Y         [-]         mass fractions
    3 : rho       [kg/m^3]    mass densities
    4 : conc      [mol/m^3]   molar concentrations
    5 : Cp        [J/mol-K]   specific heats at constant pressure
    6 : H         [J/mol]     enthalpies
    7 : S         [J/mol-K]   entropies
    8 : G         [J/mol]     Gibbs free energies
    9 : Cp        [J/kg-K]    specific heats at constant pressure
    10: H         [J/kg]      enthalpies
    11: S         [J/kg-K]    entropies
    12: G         [J/kg]      Gibbs free energies
    13: J         [kg/m^2-s]  Species diffusion fluxes (SM Ramshaw)
    14: omega     [kg/m^3-s]  production rates due to reactions
    15: Omega11   [m^2]       (1,1) pure species collision integrals
    16: Omega22   [m^2]       (2,2) pure species collision integrals
    17: Chi^h     [-]         heavy thermal diffusion ratios
    18: Dm        [m^2/s]     mixture averaged diffusion coefficients

Reaction values (same format as mixture values):
    0 : kf        [mol,m,s,K] forward reaction rate coefficients
    1 : kb        [mol,m,s,K] backward reaction rate coefficients

Other values (same format as mixture values):
    0 : Dij       [m^2/s]     multicomponent diffusion coefficients
    1 : pi_i      [-]         element potentials
    2 : N_p       [mol]       phase moles
    3 : iters     [-]         number of continuation step iterations
    4 : newts     [-]         total number of newton iterations
    5 : Fp_k      [kg/m-Pa-s] elemental diffusion fluxes per pressure gradient
    6 : Ft_k      [kg/m-K-s]  elemental diffusion fluxes per temperature gradient
    7 : Fz_k      [kg/m-s]    elemental diffusion fluxes per element mole fraction gradient
    8 : sigmaB    [S/m]       anisotropic electric conductivity
    9 : lamB_e    [W/m-K]     anisotropic electron thermal conductivity

Example:
    mppequil -T 300:100:15000 -P 101325 -m 1-3,8 air11
```
