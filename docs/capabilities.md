<a id="top"></a>

# Why Mutation++?

Glad you asked!  Mutation++ is an open-source software framework providing data, models, and algorithms to evaluate state-of-the-art thermochemical and transport properties for partially ionized gases.  The package differs from other libraries providing similar functionalities, such as [CEA](https://www.grc.nasa.gov/www/CEAWeb/), [Chemkin](https://www.ansys.com/products/fluids/ansys-chemkin-enterprise), [Cantera](https://cantera.org/), or [Kappa](https://github.com/lkampoli/kappa), in several important ways:

1. It's designed from the ground up with the full range of thermochemical nonequilibrium timescales in mind.  Equilibrium, reacting, mult-temperature, and state-to-state models can all be used seamlessly through a unique object oriented design.
2. The framework is highly extensible, allowing new models, data, or algorithms to be implemented quickly.
3. The API makes coupling to new or existing computational fluid dynamics codes simple.  C++ and Fortran interfaces are currently available, with a Python interface on the way.

## Main Features

- Thermodynamic models including the NASA 7- and 9-coefficient polynomials and the rigid-rotator / harmonic-oscillator models
- Transport coefficient models including Wilke, Gupta-Yos, and Chapmann-Enskog
- Robust, constrained multiphase equilibrium solver
- Finite-rate chemistry, including thermal nonequilibrium effects and automatic reaction type identification
- Energy exchange processes
- Gas-surface interaction mechanism
- Extensible and efficient OOP design
- Built-in command line tools to help with prepare data files and compute equilibrium mixture properties

## Still not sure?

Check out [who else is using Mutation++ in their own tools](users.md#top) or jump into our [getting started](getting-started.md#top) tutorial to get a taste of what Mutation++ can do.