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

#include "mutation++.h"
#include "TestMacros.h"
#include <catch.hpp>

using namespace Mutation;
using namespace Catch;


TEST_CASE(
    "Bug #94 fix: 1 species equilibrium cp doesn't NaN",
    "[bugs][thermodynamics]"
)
{
    MixtureOptions opts;
    opts.setSpeciesDescriptor("NH3");
    opts.setStateModel("Equil");
    Mixture mix(opts);                       

    double T = 1548;
    double P = 1106507.818257603;
    mix.setState(&P, &P, 1);

    double cp_mole = mix.mixtureEquilibriumCpMole();
    CHECK(!std::isnan(cp_mole));
    
    double cp_mass = mix.mixtureEquilibriumCpMass();
    CHECK(!std::isnan(cp_mass));

    double cp_NH3;
    mix.speciesCpOverR(&cp_NH3);
    CHECK(cp_mole == cp_NH3 * RU);
    CHECK(cp_mass == cp_NH3 * RU / mix.speciesMw(0));
}