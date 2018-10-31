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
    * VT transfer (Milikan and White (MW), Park correction)
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


