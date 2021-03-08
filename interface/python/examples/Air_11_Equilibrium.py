#!/usr/bin/env python

import mutationpp as mpp
import numpy as np
import matplotlib.pyplot as plt

from labellines import labelLines
from cycler import cycler

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

if __name__ == "__main__":
    myMixtureOptions = mpp.MixtureOptions("air_11")
    print("Mixture path:", myMixtureOptions.getSource())
    print("Mixture species:", myMixtureOptions.getSpeciesDescriptor())
    myMixtureOptions.setStateModel("Equil")
    mix = mpp.Mixture(myMixtureOptions)
    #mix = mpp.Mixture("air_11")

    ne = mix.nElements()
    ns = mix.nSpecies()

    Tin = 300.0
    Tout = 15000.0
    Temperatures = np.linspace(Tin, Tout, 1000)
    Pressure = 1000.0

    species_descriptor = np.empty((0, ns), float)
    for Temperature in Temperatures:
        mix.setState(Pressure, Temperature, 1)
        species_descriptor = np.vstack((species_descriptor, np.array(mix.X())))

    print("Average diffusion coefficients:\n", mix.averageDiffusionCoeffs())
    print("Number of elements:", mix.nElements())
    print("Number of energy equations:\n", mix.nEnergyEqns())

    print("Species molecular weight:\n", mix.speciesMw())
    print("Element matrix:\n", mix.elementMatrix())



    plt.figure()
    for i in range(ns):
        plt.plot(Temperatures, species_descriptor[:, i], label=mix.speciesName(i))
    ax = plt.gca()
    ax.yaxis.set_label_coords(0.0, 1.00)
    plt.ylabel('mole fraction', rotation=0)
    plt.xlabel('Temperature, K')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    xvals = [13000, 13000, 13000, 2500, 400, 5500, 6000, 7000, 8000, 400, 400]
    labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left', xvals=xvals)
    plt.savefig('composition.eps', transparent=True)

    mix.addComposition("N2:1.0, O2:0.0", True)

    species_descriptor = np.empty((0, ns), float)
    for Temperature in Temperatures:
        mix.setState(Pressure, Temperature, 1)
        species_descriptor = np.vstack((species_descriptor, np.array(mix.X())))

    plt.figure()
    for i in range(ns):
        plt.plot(Temperatures, species_descriptor[:, i], label=mix.speciesName(i))
    ax = plt.gca()
    ax.yaxis.set_label_coords(0.0, 1.00)
    plt.ylabel('mole fraction', rotation=0)
    plt.xlabel('Temperature, K')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    xvals = [13000, 13000, 13000, 2500, 400, 5500, 6000, 7000, 8000, 400, 400]
    labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left', xvals=xvals)
    plt.savefig('composition_n2.eps', transparent=True)

    plt.show()
