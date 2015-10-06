/**
 * @file MixtureFixtures.h
 *
 * @brief Provides test fixtures which provide preloaded mixtures for different
 * test-cases to use.
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

#ifndef TESTS_MIXTURE_FIXTURES_H
#define TESTS_MIXTURE_FIXTURES_H

#include <boost/test/included/unit_test.hpp>
#include "mutation++.h"

class MixtureFromCommandLine {
public:
    MixtureFromCommandLine()
        : mp_mix(NULL)
    {
        BOOST_REQUIRE_MESSAGE(
            boost::unit_test::framework::master_test_suite().argc > 1,
            "Mixture name must be given in command line arguments");

        // Load the mixture from the first command line argument
        BOOST_TEST_MESSAGE("Loading mixture from command line argument");
        mp_mix = new Mutation::Mixture(
            boost::unit_test::framework::master_test_suite().argv[1]);
    }

    ~MixtureFromCommandLine()
    {
        // Load the mixture from the first command line argument
        BOOST_TEST_MESSAGE("Destroying mixture from command line argument");
        delete mp_mix;
    }

    Mutation::Mixture& mix() { return *mp_mix; }

private:
    Mutation::Mixture* mp_mix;
};

#endif // TESTS_MIXTURE_FIXTURES_H
