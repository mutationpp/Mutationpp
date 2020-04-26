<a id="top"></a>

# Release notes
**Contents**<br>
[1.0.1](#101)<br>
[1.0.0](#100)<br>
[Older versions](#older-versions)<br>

## 1.0.1

**Improvements**
- Added NH3 and CH2O transport data
- Added heavy thermal diffusion ratios to fortran wrapper

**Fixes**
- Fixed Parsing of STS ionized species
- Fixed bugs and segfault in species chemistry jacobian calculation
- Fixed convergence error EquilStateModel::setState() for mixture energies near zero
- Fixed FPE or NaN occuring with computing equilibrium Cp for a single species mixture
- Fixed missing include causing compiler errors on Mac OSX
- Fixed wrong type name used Transport.cpp which was silently failing
- Fixed segfault when computing pure species collision integrals for mixtures without electrons 

**Miscellaneous**
- Updated copyright year
- Switched from Travis to Github Actions for CI
- Documentation updates

## 1.0.0

This is the first "official" version of Mutation++.  All future releases should maintain these release notes by listing the `Improvements`, `Fixes`, and `Miscellaneous` updates from the previous version.

## Older versions

Older versions of the library exist but unfortunately did not maintain any release notes.
