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
using namespace Mutation::Transfer;
using namespace Catch;
using namespace Eigen;


TEST_CASE("Energy transfer source terms are zero in equilibrium",
    "[transfer]"
)
{
    const std::string transfer_terms [] = {
        "OmegaCE", "OmegaCElec", "OmegaCV", "OmegaET", "OmegaI", "OmegaVT"
    };

    const double tol = std::numeric_limits<double>::epsilon()*1000;

    MIXTURE_LOOP
    (
        // Loop over the different transfer terms
        for (int i = 0; i < 6; ++i) {
            SECTION(transfer_terms[i]) {
                TransferModel* p_omega =
                     Utilities::Config::Factory<TransferModel>::create(
                         transfer_terms[i], mix);

                EQUILIBRATE_LOOP
                (
                    double rhoe = mix.mixtureEnergyMass()*mix.density();
                    CHECK(p_omega->source() / rhoe == Approx(0.0).margin(tol));
                )

                delete p_omega;
            }
        }

        SECTION("Total") {
            VectorXd omega(mix.nEnergyEqns()-1);

            EQUILIBRATE_LOOP
            (
                mix.energyTransferSource(omega.data());
                double rhoe = mix.mixtureEnergyMass()*mix.density();
                for (int i = 0; i < mix.nEnergyEqns()-1; ++i)
                    CHECK(omega[i]/rhoe == Approx(0.0).margin(tol));
            )
        }
    )
}

