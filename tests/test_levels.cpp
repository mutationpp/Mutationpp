/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include <catch.hpp>
#include <Eigen/Dense>

using namespace Mutation;
using namespace Catch;
using namespace Eigen;
using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

void testOnce(std::string specName)
{
    // setup the mixture
    MixtureOptions mixture_options;
    mixture_options.setThermodynamicDatabase("RRHO");
    mixture_options.setSpeciesDescriptor(specName+"(*)");

    Mixture mix(mixture_options);

    const int ns = mix.nSpecies();

    Eigen::ArrayXd xi(ns);
    Eigen::ArrayXd xth(ns);
    Eigen::ArrayXd degen(ns);
    Eigen::ArrayXd en(ns);

    // get the energy levels
    IO::XmlDocument species_doc(databaseFileName("species.xml", "thermo"));

    ParticleRRHO spec(*(species_doc.root().findTagWithAttribute(
        "species", "name", specName)->
            findTagWithAttribute("thermodynamics", "type", "RRHO")));

    CHECK(spec.nElectronicLevels() == ns);

    for (int i = 0; i < ns; ++i) {
        en[i] = spec.electronicEnergy(i).second;
        degen[i] = spec.electronicEnergy(i).first;
    }

    // check that levels follow the Boltzmann distribution at equilibrium
    EQUILIBRATE_LOOP
    (
        // get the species mole fractions
        xi = Eigen::Map<const Eigen::ArrayXd>(mix.X(), ns);

        double t1 = -1.0/T;

        // compute Boltzmann mole fractions
        xth = degen*( (en*t1).exp() );
        xth = xth/(xth.sum());
        CHECK(xth.sum() == Approx(1.0));

        for (int i = 0; i < ns; ++i)
            CHECK( xi[i] == Approx(xth[i]) );
    )

}

TEST_CASE("Electronic levels are loaded properly using star token",
    "[thermodynamics]"
)
{
    testOnce("Ar");
    testOnce("N");
    testOnce("O2");
}


