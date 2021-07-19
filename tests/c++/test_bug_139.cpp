/*
 * Copyright 2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <memory>

#include "mutation++.h"
#include <Eigen/Dense>
#include <catch.hpp>

using namespace Mutation;
using namespace Catch;
using namespace Eigen;

/*

This is the test data provided in the issue:

air_5 - RRHO 
rhos: {0.0831209815439599, 0.03252762384728387, 0.0005480571702340367, 0.02473956358534552, 8.688041095193619e-06}
input T: 8042.055877561177
input Tv: 14387.97157416015
input rho_e: 4877971.072365317
input rho_eV: 673333.7559139833
expected e: 34609060.57143456
returned e: 34609060.57143456
expected eV: 4777283.095268293
returned eV: 4777283.095268293
returned T: 8042.055877561177
returned Tv: 14387.97157416015
 
air_11 - RRHO 
rhos: {3.772083171225316e-07, 0.008217402454188135, 0.00158203164044893, 3.349905442259246e-06, 3.818030661535437e-06, 5.411162803271429e-08, 0.03034510868639386, 0.01013019524997852, 1.718604852728191e-06, 1.311304142186017e-05, 6.967273753567727e-08}
input T: 11742.36111897272
input Tv: 5078.056874545492
input rho_e: 2976361.926591469
input rho_eV: 80226.49707629484
expected e: 59175453.94295745
returned e: 59175453.94297452
expected eV: 1595047.746152252
returned eV: 2336397.410620507
returned T: 11076.42304153593
returned Tv: 6539.047348525341

 
CO2_8 - RRHO 
rhos: {1.634350644243097e-07, 0.003106515000745763, 0.0006160141919509619, 0.02351932703393298, 0.07072164609821735, 0.0007229419681532173, 3.470025039233278e-08, 1.019784785257949e-05}
input T: 15212.22216121357
input Tv: 9718.179857642906
input rho_e: 4375416.313420466
input rho_eV: 138236.0407565047
expected e: 44331878.31725346
returned e: 44331878.31725392
expected eV: 1400612.627209754
returned eV: 1672067.877769741
returned T: 14925.69634400333
returned Tv: 10843.8814584733

*/

// Just create a mixture for our purposes
std::shared_ptr<Mixture> getMixture(const std::string& species)
{
    MixtureOptions opts;
    opts.setSpeciesDescriptor(species);
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    return std::shared_ptr<Mixture>(new Mixture(opts));
};


// Test a single mixture
void testMixture(
    const std::string& species, const std::vector<double>& rhos,
    const std::vector<double>& Ts)
{
    auto mix = getMixture(species);
    mix->setState(rhos.data(), Ts.data(), 1);
    double rho = mix->density();
    
    ArrayXd input_rhoe(mix->nEnergyEqns());
    mix->mixtureEnergies(input_rhoe.data());
    input_rhoe *= rho;

    mix->setState(rhos.data(), input_rhoe.data(), 0);
    
    ArrayXd output_rhoe(mix->nEnergyEqns());
    mix->mixtureEnergies(output_rhoe.data());
    output_rhoe *= rho;

    for (int i = 0; i < mix->nEnergyEqns(); ++i) {
        std::cout << input_rhoe[i] << " == " << output_rhoe[i] << '\n';
        CHECK(input_rhoe[i] == output_rhoe[i]);
    }
}


// Check the three cases associated with issue 139
TEST_CASE(
    "Bug #139 fix: vibronic energies unchanged after setState()",
    "[bugs][thermodynamics]"
)
{
    testMixture("N O NO N2 O2", {
            0.0831209815439599, 0.03252762384728387, 0.0005480571702340367, 
            0.02473956358534552, 8.688041095193619e-06
        }, {
            8042.055877561177, 14387.97157416015
        }
    );

    testMixture("e- N+ O+ NO+ N2+ O2+ N O NO N2 O2", {
            3.772083171225316e-07, 0.008217402454188135, 0.00158203164044893, 
            3.349905442259246e-06, 3.818030661535437e-06, 5.411162803271429e-08,
            0.03034510868639386, 0.01013019524997852, 1.718604852728191e-06, 
            1.311304142186017e-05, 6.967273753567727e-08
        }, {
            11742.36111897272, 5078.056874545492
        }
    );

    testMixture("e- C+ O+ C O CO CO2 O2", {
            1.634350644243097e-07, 0.003106515000745763, 0.0006160141919509619, 
            0.02351932703393298, 0.07072164609821735, 0.0007229419681532173, 
            3.470025039233278e-08, 1.019784785257949e-05
        }, {
            15212.22216121357, 9718.179857642906
        }
    );
}