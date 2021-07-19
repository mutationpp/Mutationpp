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

O2 Dissociation with the Mutation++ python biding

```python
# Import desired packages
import mutationpp as mpp
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines
from cycler import cycler
```

Load the mutation++ package

```python
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

```python
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

```python
opts = mpp.MixtureOptions("Mars_19")
opts.setThermodynamicDatabase("NASA-9")
opts.setStateModel("ChemNonEq1T")
print("Mixture species:", opts.getSpeciesDescriptor())
```

```python
mix = mpp.Mixture(opts)
ns = mix.nSpecies()
```

```python
T = 300
P = 10000
```

```python
mix.equilibrate(T,P)
```

```python
rhoi_eq = mix.densities()
mix.setState(rhoi_eq,T,1)
```

```python
Tinit = 15000
mix.setState(rhoi_eq,Tinit,1)
total_energy = mix.mixtureEnergyMass()*mix.density()
```

```python
tol = 1.e-7
conv = 1.
time = 0.
time_final = 1e-7
dt = 1.e-9
rhoi = rhoi_eq
```

```python
y0 = np.divide(rhoi_eq,mix.density())
x0 = mix.convert_y_to_x(y0)
temperature = np.array(mix.T())
mass_fractions = np.array([y0])
mole_fraction =  np.array([x0])
time_reaction = np.array(time)
```

```python
while (time < time_final):
    dt = min(dt * 1.0002, 1.e-10)
    drhoi = runge_kutta4(mix, rhoi, total_energy, dt)
    rhoi += drhoi
    conv = max(abs(np.array(drhoi)))
    time += dt
    mix.setState(rhoi,total_energy,0)
    y = mix.Y()
    x = mix.convert_y_to_x(np.array(y))
    mass_fractions = np.vstack((mass_fractions,[y]))
    mole_fraction =  np.vstack((mole_fraction,[x]))
    temperature = np.vstack((temperature,mix.T()))
    time_reaction  = np.vstack((time_reaction,time))
```

```python
plt.figure()
for i in range(ns):
    plt.plot(time_reaction,mole_fraction[:,i],label=mix.speciesName(i))
ax = plt.gca()
ax.yaxis.set_label_coords(0.0, 1.00)
plt.ylabel('mass fraction', rotation=0)
plt.xlabel('time, s')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)  
labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left')
```

```python
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

```python

```
