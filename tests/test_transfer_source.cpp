/**
 * @file test_transfer_source.cpp
 *
 * @brief Contains tests related to energy transfer source terms.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE transfer source tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <numeric>
#include "MixtureFixtures.h"
#include "AutoRegistration.h"
#include "TransferModel.h"

using namespace Mutation;
using namespace Mutation::Transfer;
using namespace Mutation::Utilities::Config;

// Provide mixture fixture to all test cases in this suite
BOOST_FIXTURE_TEST_SUITE(Omega, MixtureFromCommandLine)

/**
 * Small helper function to compute equilibrium mixtures and make sure the given
 * energy transfer source is zero at several temperatures and pressures.
 */
void checkOmegaIsZeroInEquil(
    MixtureFromCommandLine* p_fix, const std::string& name)
{
    Mixture& mix = p_fix->mix();
    TransferModel* p_omega = Factory<TransferModel>::create(name, mix);
    const double tol = 1.0e-8 * mix.nSpecies();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state
            p_fix->equilibrateMixture(T, P);

            // Make sure the source term is near zero
            double rhoe = mix.mixtureEnergyMass()*mix.density();
            BOOST_CHECK_SMALL(p_omega->source()/rhoe, tol);
        }
    }

    // Clean up the source term from memory
    delete p_omega;
}

BOOST_AUTO_TEST_CASE(OmegaCE_is_zero_in_equil)
{
    if (mix().hasElectrons())
    	checkOmegaIsZeroInEquil(this, "OmegaCE");
}

BOOST_AUTO_TEST_CASE(OmegaCElec_is_zero_in_equil)
{
	checkOmegaIsZeroInEquil(this, "OmegaCElec");
}

BOOST_AUTO_TEST_CASE(OmegaCV_is_zero_in_equil)
{
    checkOmegaIsZeroInEquil(this, "OmegaCV");
}

BOOST_AUTO_TEST_CASE(OmegaET_is_zero_in_equil)
{
    if (mix().hasElectrons())
        checkOmegaIsZeroInEquil(this, "OmegaET");
}


BOOST_AUTO_TEST_CASE(OmegaI_is_zero_in_equil)
{
	if (mix().hasElectrons())
	    checkOmegaIsZeroInEquil(this, "OmegaI");
}

BOOST_AUTO_TEST_CASE(OmegaVT_is_zero_in_equil)
{
    checkOmegaIsZeroInEquil(this, "OmegaVT");
}

/*
 * Ensures that the total energy transfer source terms associated with the
 * StateModel is zero when mixture is in thermochemical equilibrium for a range
 * of temperatures and pressures.
 */
BOOST_AUTO_TEST_CASE(Omega_is_zero_in_equil)
{
    std::vector<double> omega(mix().nEnergyEqns(), 0.0);
    std::vector<double> em(mix().nEnergyEqns(), 0.0);

    // Tolerance must account for errors in equilibrium solution which are
    // propagated forward into the production rates
    const double tol = 1.0e-8 * mix().nSpecies();

    // Loop over range of pressures and temperatures
    for (int ip = 0.0; ip < 10; ++ip) {
        double P = std::exp(ip/9.0*std::log(100000.0)+std::log(10.0));
        for (int it = 0.0; it < 10; ++it) {
            double T = 1000.0*it + 1000.0;

            // Set an equilibrium state
            equilibrateMixture(T, P);
            double rhoe = mix().density()*mix().mixtureEnergyMass();

            // Make sure the energy transfer sources are close to zero
            mix().energyTransferSource(&omega[0]);
            for (int i = 0; i < mix().nEnergyEqns()-1; ++i)
                BOOST_CHECK_SMALL(omega[i]/rhoe, tol);
        }
    }
}

// End of the transfer source test suite
BOOST_AUTO_TEST_SUITE_END()
