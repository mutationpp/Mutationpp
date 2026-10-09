# Changelog
All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project adheres to [Semantic Versioning](http://semver.org/). For each version listed below, changes are grouped to describe their impact on the project, as follows:

- `Added` for new features
- `Changed` for changes in existing functionality
- `Deprecated` for once-stable features removed in upcoming releases
- `Removed` for deprecated features removed in this release
- `Fixed` for any bug fixes
- `Security` to invite users to upgrade in case of vulnerabilities

## [Unreleased]

## [1.3.0] - 2026-07-23
### Added  
 * `Energy exchange mechanisms`
    * VT transfer (Millikan and White (MW), Park correction)
        - Enabling to select original mixing rule or Gnoffo at runtime     
 * `GSI`
    * Functions to get surface emissivity and radiative heat flux
    * Function to set GSI solver convergence tolerance
 * `Examples`
    * 2T dissociation of O2

## [1.2.0] - 2026-06-15
### Changed
 * `Included Tools`
    * bprime
        - Enabling to pass surface composition
  
## [1.1.3] - 2023-08-23
### Fixed
 * `Thermodynamics`:
    * Pure species thermo: 
        * Update electron molecular weight

## [1.1.2] - 2022-01-06
### Fixed
 * `Interfaces to other languages`
    * Fix issue with Fortran compiler flags not being set
        
## [1.1.1] - 2021-09-24
### Fixed
 * `Interfaces to other languages`
    * Fix install path for fortran module and export name
        
## [1.0.5] - 2021-07-19
### Added
 * `Interfaces to other languages`
    * Python

## [1.0.1] - 2020-04-26
### Added
 * `Thermodynamics`:
    * Databases: 
        * NASA-7, 
        * NASA-9, 
        * RRHO (custom for some species)
    * Pure species thermo
    * Mixture thermo (mostly multi-temperature models)

* `Transport`
    * Custom collision integral database which supports:
        * Tables
        * Curve-fits from the Capitelli group
        * Old Mutation-style curve-fits (NASA format)
        * Constants, ratios of other integrals
    * Multiple algorithms for diffusion, thermal conductivity, thermal diffusion ratios, viscosity, electric conductivity 
    * Magnetized transport
  
* `Gas phase chemistry`
    * Mechanism sanity check
    * Arrhenius rate laws
    * Elementary reactions
    * Third-body reaction special treatment
    * Multi-temperature rate laws
  
* `GSI`
    * Surface Mass Balance
        * Gamma models
    * Surface Energy Balance
  
* `Energy exchange mechanisms`
    * VT transfer (Millikan and White (MW), Park correction)
    * ET transfer
    * Chemical energy exchange:
        * Electron impact ionization and dissociation
        * Vibration-chemistry-vibration (simple non-preferential dissociation)
        * Chemistry-electronic energy coupling

* `State Models`
    * Equil: thermochemical equilibrium
    * ChemNonEq1T: thermal equilibrium, chemical nonequilibrium
    * ChemNonEqTTv: two temperature chemical nonequilibrium (Park model)

* `Included Mixtures/Models/Data`
  
    * Air5
    * Air11
    * Air13
    * CO28
   
* `Interfaces to other languages`
    * Fortran

* `Supported OSs`
    * Linux
    * Mac OS X

* `Included Tools`
    * checkmix
    * mppequil
    * bprime
    * mppcalc
    * mppshock


