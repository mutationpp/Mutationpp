---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region slideshow={"slide_type": "slide"} -->
$$\require{mhchem}$$   

# Adiabatic reactor for an 11 species air mixture

In this notebook, we will use mutationpp package to compute an adiabatic reactor for an air mixture composed by eleven species such that,
$$
\mathcal{S} = \{\ce{e-}, \ce{ N+},\ce{ O+}, \ce{ NO+}, \ce{ N2+} , \ce{ O2+}, \ce{ N}, \ce{ O}, \ce{ NO}, \ce{ N2}, \ce{ O2} \}.
$$

For this case, the species continuity equation reduces to an ODE and the mixture energy and pressure are constant,

$$
\partial_t \rho_i = \dot{\omega_i} \quad \forall i \in \mathcal{S}, \\
\rho e(\rho_i, T) = \text{const}.
$$

We have integrated the ODE with the Runge–Kutta 4th (**RK4**) order method (https://en.wikipedia.org/wiki/Runge–Kutta_methods)

The mixture energy is initialized at $T_{\text{init}}$ = 15000 K and P = 10000 Pa, and the initial composition is obtained at the equilibrium conditions of $T_{\text{eq}}$ = 300 K and P = 10000 Pa.

*Run the following cell to initialize the necessary packages*
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
# Import desired packages
import mutationpp as mpp
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines
from cycler import cycler
```

<!-- #region slideshow={"slide_type": "slide"} -->
In the following cell we define the  Runge–Kutta function used to integrate the species continuity equation.
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
def runge_kutta4(mix, rhoi, T, dt):
    
    mix.setState(rhoi,T,0)
    wdot_1 = mix.netProductionRates() 

    rhoi += 0.5 * np.multiply(dt, wdot_1)
    mix.setState(rhoi,T,0)
    wdot_2 = mix.netProductionRates() 
    
    rhoi += 0.5 * np.multiply(dt, wdot_2)
    mix.setState(rhoi,T,0)
    wdot_3 = mix.netProductionRates() 
    
    rhoi += np.multiply(dt,wdot_3)
    mix.setState(rhoi,T,0)
    wdot_4 = mix.netProductionRates() 

    return 1./6. * np.multiply(dt, (np.array(wdot_1) + 2 * np.array(wdot_2) + 2 * np.array(wdot_3) + np.array(wdot_4)))
```

<!-- #region slideshow={"slide_type": "notes"} -->
In this cell, we define some parameters for beautiful plotting.
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 18

plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
cc = cycler('linestyle', ['-', '--', ':', '-.', '-']) * cycler('color', ['r', 'g', 'b', 'm', 'k'])
plt.rc('axes', prop_cycle=cc)

plt.rc('text', usetex=True)
plt.rc('font', family='Helvetica')
```

<!-- #region slideshow={"slide_type": "slide"} -->
## Initiallizing the mixture from mutation++

The magic starts in the following cells. 

Before initializing the mixture, we need to set-up the desired options for our adiabatic conditions.
The first step is to specify the mixture name, which is the **air_11** mixture, and the thermal database we want to use (**NASA-9**).
Afterward, we specify the state model to be used on our reactor. In this case, we are interested in a one-temperature chemical non-equilibrium model. 
Other models can be used, such as **ChemNonEqTTv** for a two-temperature model and **EquilTP** for a thermo-chemical equilibrium model.

We can print the species within our mixture by using the **getSpeciesDescriptor()** method.
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
opts = mpp.MixtureOptions("air_11")
opts.setThermodynamicDatabase("NASA-9")
opts.setStateModel("ChemNonEq1T")
print("Mixture species:", opts.getSpeciesDescriptor())
```

<!-- #region slideshow={"slide_type": "slide"} -->
After, we initialize our mixture based on the options chosen above.

Note, that one could initilize the mixture directly using the mixture name such that,

```python
     mix = mpp.Mixture("air_11")
```

The options above would need to be specified directly in the mixture file **air_11.xml**.
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
mix = mpp.Mixture(opts)
ns = mix.nSpecies()
```

<!-- #region slideshow={"slide_type": "slide"} -->
We define the initial conditions of our mixture at $T_{\text{eq}}$ = 300 K and P = 10000 Pa using the method equilibrate.
The partial densities are obtained with the method mix.densities().
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
T = 300
P = 10000
mix.equilibrate(T,P)
rhoi_eq = mix.densities()
```

<!-- #region slideshow={"slide_type": "slide"} -->
Later, we set the state using the intial composition and temperature $T_{\text{init}}$ = 15000.


*Note that one could set state using different options such that:*
  * Sets the current state of the mixture.  Variable sets can be the following:
     *   0: conserved variables (species densities, total energy density)
     *   1: primitive set 1 (species densities, mixture temperature)
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
Tinit = 15000
mix.setState(rhoi_eq,Tinit,1)
total_energy = mix.mixtureEnergyMass()*mix.density()

```

<!-- #region slideshow={"slide_type": "slide"} -->
We now define the initial and finaly conditions of our RK4 method:
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
time = 0.
time_final = 1e-7
dt = 1.e-9
rhoi = rhoi_eq

y0 = np.divide(rhoi_eq,mix.density())
x0 = mix.convert_y_to_x(y0)
temperature = np.array(mix.T())
mass_fractions = np.array([y0])
mole_fraction =  np.array([x0])
time_reaction = np.array(time)
```

<!-- #region slideshow={"slide_type": "slide"} -->
## ODE integration 

After everything is set, we can start integrating our ODE.
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
while (time < time_final):
    dt = min(dt * 1.0002, 1.e-10)
    drhoi = runge_kutta4(mix, rhoi, total_energy, dt)
    rhoi += drhoi
    time += dt
    mix.setState(rhoi,total_energy,0)
    y = mix.Y()
    x = mix.convert_y_to_x(np.array(y))
    mass_fractions = np.vstack((mass_fractions,[y]))
    mole_fraction =  np.vstack((mole_fraction,[x]))
    temperature = np.vstack((temperature,mix.T()))
    time_reaction  = np.vstack((time_reaction,time))
```

<!-- #region slideshow={"slide_type": "slide"} -->
## Plotting the solution

Finally, we can plot the evolution of the species mole fraction and the temperature relaxation in time
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
plt.figure()
for i in range(ns):
    plt.plot(time_reaction,mole_fraction[:,i],label=mix.speciesName(i))
ax = plt.gca()
ax.yaxis.set_label_coords(0.0, 1.00)
plt.ylabel('mole fraction', rotation=0)
plt.xlabel('time, s')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)  
labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left')
```

```python slideshow={"slide_type": "subslide"}
plt.figure()
plt.plot(time_reaction,temperature)
ax = plt.gca()
ax.yaxis.set_label_coords(0.0, 1.00)
plt.ylabel('Temperature, K', rotation=0)
plt.xlabel('time, s')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)  
labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left')
```
