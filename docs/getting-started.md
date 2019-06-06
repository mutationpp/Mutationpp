<a id="top"></a>

# Oops! We have some work left to do on this page...

**Outline** <br>
- Installation, point to installation instructions
- Getting started with simple mixture
  - Loading mixture
  - Printing some mixture values
  - Equilibrating the mixture
  - Setting the state
  - Printing some mixture properties
- Changing mixture options with MixtureOptions object
- Building your own mixture models, point to input file docs
- Command line tools
- Going further, provided examples, reference section, etc


## Setting up

If you want to follow along with the tutorial below, be sure to work through the [installation instructions](installation.md#top) before proceeding.  If you just want a quick look at how Mutation++ works, that's cool too!

## Getting started with a simple mixture

Mutation++ exposes most of its functionality through a single class called `Mixture`.  A `Mixture` object represents an entire thermochemical model for a gas mixture.

### Loading a mixture

Let's create a new mixture object, using a simple 5 species air model.

```c++
Mixture mix("air_5");
```

Note that mixture models can be referred to by name.  In this case, we used "air_5".  By default, several mixtures are included with Mutation++.  We will talk more about how to create your own mixtures in a little bit, but for now let's continue with this one.

### Print some properties of the mixture model

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

