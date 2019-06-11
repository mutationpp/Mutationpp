<a id="top"></a>

# Getting Started

This simple tutorial is a good way to start coding with Mutation++ and will walk you through the basics of the library.

## Setting up

If you want to follow along with the tutorial below, be sure to work through the [installation instructions](installation.md#top) before proceeding.  The full code of this tutorial can be found in the `examples/c++/tutorial` directory of the project.

If you just want a quick look at how Mutation++ works, that's cool too!

## Getting started with a simple mixture

Mutation++ exposes most of its functionality through a single class called `Mixture`.  A `Mixture` object represents an entire thermochemical model for a gas mixture.

### Loading a mixture

You can instantiate a new mixture by simply using a defined mixture name.

```c++
Mixture mix("air_5");
```

In this case, we used "air_5".  By default, several mixtures are included with Mutation++.  We will talk more about how to create your own mixtures in a little bit, but for now let's continue with this one.

### Check model parameters

Now that we have our mixture loaded, let's print a few statistics about the model.

**Code**
```c++
std::cout << "# of elements: " << mix.nElements() << '\n';
std::cout << "# of species: " << mix.nSpecies() << ' ';
std::cout << mix.nGas() << " (gas) " << mix.nCondensed() << " (condensed)\n";
std::cout << "# of reactions: " << mix.nReactions() << '\n';
std::cout << "# of temperatures: " << mix.nEnergyEqns() << '\n';
```

**Output**
```
# of elements: 2
# of species: 5 5 (gas) 0 (condensed)
# of reactions: 5
# of temperatures: 1
```

We see that the air_5 model has 5 gaseous species comprised of two elements and 5 reactions.  We can go a little deeper and print the species names and reaction formulas as well.

**Code**
```c++
std::cout << "Species:\n";
for (auto& s: mix.species()) 
    std::cout << s.name() << ' ';
std::cout << "\n\n";

std::cout << "Reactions:\n";
for (auto& r: mix.reactions())
    std::cout << r.formula() << '\n';
std::cout << '\n';
```

**Output**
```
Species:
N O NO N2 O2 

Reactions:
N2+M=2N+M
O2+M=2O+M
NO+M=N+O+M
N2+O=NO+N
NO+O=O2+N
```

### Setting the mixture state

The primary purpose of the `Mixture` class is to allow the user to alter the thermodynamic state of the mixture, and to compute various properties of the mixture at the given state.  The thermodynamic state will depend on the type of `StateModel` used by the mixture, but more on that later.  In general, we can change the state of the mixture in two ways.

1. Equilibrate the mixture at a given temperature and pressure, or
2. Set the thermodynamic variables explicitly.

Equilibrating the mixture can be done easily with a call to `equilibrate`.

```c++
mix.equilibrate(300.0, ONEATM);
```

The above sets the mixture composition to the equilibrium composition at with a fixed temperature of 300 K and pressure of 101,325 Pa or 1 atm.  Default elemental mole fractions are also used, but this could be changed with an additional parameter.  Note the use of the defined constant for one atmosphere.

Let's check the mixture state after this call by printing the temperature, pressure, and species mole fractions.

**Code**
```c++
std::vector<double> T(mix.nEnergyEqns());
mix.getTemperatures(T.data());

std::cout << "\nTemperature(s) [K]: ";
for (auto v: T) std::cout << v << ' ';

std::cout << "\nPressure [Pa]: " << mix.P();

std::cout << "\nSpecies mole fractions:\n";
for (int i = 0; i < mix.nSpecies(); ++i)
    std::cout << mix.speciesName(i) << ": " << mix.X()[i] << '\n';
```

**Output**
```
Temperature(s) [K]: 300 
Pressure [Pa]: 101325
Species mole fractions:
N: 4.94091e-80
O: 2.26658e-41
NO: 2.59341e-16
N2: 0.79
O2: 0.21
```

Note that it is possible for a mixture to have multiple temperatures, depending on the state model.  We have allowed for this in the code above, even though, in this example, there is only a single temperature.

Equilibrating the mixture at a higher temperature,

```c++
mix.equilibrate(5000.0, ONEATM);
```

yields a different composition.

```
Temperature(s) [K]: 5000 
Pressure [Pa]: 101325
Species mole fractions:
N: 0.0269021
O: 0.324738
NO: 0.0174998
N2: 0.628901
O2: 0.00195878
```

We see at this temperature, the oxygen molecules have mostly dissociated into atoms or combined with nitrogen atoms to form NO.

As previously mentioned, we can also set the state by explicitly providing the state variables.  This is done using the `setState` method of the mixture object. This function may have different behaviors, depending on the state model currently loaded.  In general however, `setState` takes 3 arguments: 

1. a "mass" array,
2. an "energy" array, and
3. an integer specifying how the state arrays should be interpreted.

Let's provide a concrete example.  First, let's get the mixture's current temperatures and mass densities from the previous equilibrate call.

```c++
std::vector<double> rho(mix.nSpecies()), T(mix.nEnergyEqns());
mix.getTemperatures(T.data());
mix.densities(rho.data());
```

Now, we let's explicitly set the state, using species densities and mixture temperatures.  For the default state model used in the `air_5` example, this corresponds to an index of 1.

```c++
const int use_rhoi_and_t = 1;
mix.setState(rho.data(), T.data(), use_rhoi_and_t);
```

Printing the mixture state as before, we see that is unchanged, as expected.

```
Temperature(s) [K]: 5000 
Pressure [Pa]: 101325
Species mole fractions:
N: 0.0269021
O: 0.324738
NO: 0.0174998
N2: 0.628901
O2: 0.00195878
```

While we used the current state to demonstrate the `setState` call here, typically the state information will be known.  For example, this could be provided by an external CFD tool which is coupled to Mutation++.

### Computing physicochemical properties

Once we have set the state of the mixture, we can proceed to compute a variety of thermochemical properties of the gas.  For example, here are a few properties computed with the state set in the previous section.

**Code**
```c++
std::cout << "\nMixture enthalpy [J/kg]: " << mix.mixtureHMass();
std::cout << "\nMixture entropy [J/kg-K]: " << mix.mixtureSMass();

std::vector<double> lambda(mix.nEnergyEqns());
mix.frozenThermalConductivityVector(lambda.data());
std::cout << "\nMixture thermal conductivities [W/m-K]: ";
for (auto v: lambda) std::cout << v << ' ';

std::cout << "\nMixture dynamic viscosity [Pa-s]: " << mix.viscosity();
```

**Output**
```
Mixture enthalpy [J/kg]: 9.9924e+06
Mixture entropy [J/kg-K]: 11335.9
Mixture thermal conductivities [W/m-K]: 0.280127 
Mixture dynamic viscosity [Pa-s]: 0.000143242
```

Of course, there are many more properties implemented in Mutation++.  See the API documentation for more details.

## Specifying non-default mixture options

In the above examples, we loaded the mixture using only the mixture name.  This name refers to a mixture file in the data directory provided with the library.  Sometimes, it is necessary to modify the default options listed in the mixture file.  This can be done using the `MixtureOptions` class. 

```c++
MixtureOptions opts("air_5");
opts.setStateModel("ChemNonEqTTv");
opts.setThermodynamicDatabase("RRHO");
opts.setViscosityAlgorithm("Gupta-Yos");
Mixture mix(opts);
```

Note that we can instantiate a mixture options object using the [default options](input-files.md#mixture-options) in a mixture file.  You can also create a MixtureOptions object with no defaults, but you must specify every option manually in this case.  In the above, we have used the defaults for the `air_5` mixture, but explicitly set the state model, thermodynamic database, and viscosity algorithm.  For example, we have replaced the 1-temperature default model with the 2-temperature nonequilibrium model.  Printing the mixture properties with the new options yields the following output.

```
# of elements: 2
# of species: 5 5 (gas) 0 (condensed)
# of reactions: 5
# of temperatures: 2

Species:
N O NO N2 O2 

Reactions:
N2+M=2N+M
O2+M=2O+M
NO+M=N+O+M
N2+O=NO+N
NO+O=O2+N

Temperature(s) [K]: 5000 5000 
Pressure [Pa]: 101325
Species mole fractions:
N: 0.0269021
O: 0.324738
NO: 0.0174998
N2: 0.628901
O2: 0.00195878

Mixture enthalpy [J/kg]: 9.9924e+06
Mixture entropy [J/kg-K]: 11335.9
Mixture thermal conductivities [W/m-K]: 0.237793 0.0423338 
Mixture dynamic viscosity [Pa-s]: 0.000143242
```

As we can see, the number of temperatures in the model is now 2.  We can also see that there are now 2 thermal conductivities.  These correspond to the translation-rotation and vibrational-electronic-electron energy modes used in the 2T model.

## Building your own mixture models

The easiest way to construct a new mixture model is to modify one of the existing models that comes packaged with Mutation++.  See the documentation on [input files](input-files.md#top) for more information.

##  Built-in command-line tools

A number of tools are provided out-of-the-box with Mutation++.  For example, the [checkmix](checkmix.md#top) tool can be used to print information about a given mixture file, and to check that all data is loaded correctly.  This is very useful if you want to see what models are currently available, or to check if a custom model is implemented correctly.  The [mppequil](mppequil.md#top) tool computes equilibrium properties for a given mixture over a prescribed temperature and pressure range.  Check the [reference section](Readme.md#top) for an overview of all the tools provided.

## Going further

That's it for this simple tutorial.  If you want to learn more, take a look at the full [reference section](Readme.md#top) for all the documentation on Mutation++.