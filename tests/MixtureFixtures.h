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

    /**
     * Sets the mixture to an equilibrium state at the given T and P.  Assumes
     * that the variable set 1 of the state model accepts species densities and
     * temperature vectors.
     */
    void equilibrateMixture(double T, double P)
    {
        static std::vector<double> rhoi(mix().nSpecies());
        static std::vector<double> tmps(mix().nEnergyEqns());

        // Get equilibrium mole fractions
        mix().equilibriumComposition(T, P, &rhoi[0]);

        // Convert to species densities
        double rho = mix().density(T, P, &rhoi[0]);
        mix().convert<Mutation::Thermodynamics::X_TO_Y>(&rhoi[0], &rhoi[0]);
        for (int i = 0; i < mix().nSpecies(); ++i)
            rhoi[i] *= rho;

        tmps.assign(mix().nEnergyEqns(), T);
        mix().setState(&rhoi[0], &tmps[0], 1);
    }

private:
    Mutation::Mixture* mp_mix;
};

/*
 *
 */

 class SetMixtureGSIO2gamma {
 public:
     SetMixtureGSIO2gamma()
         : mp_mix(NULL)
     {
         BOOST_TEST_MESSAGE("Loading mixture for gsi test");
         mp_mix = new Mutation::Mixture("O2_gsi_gamma");
 
         rhoi.resize(mix().nSpecies());
         tmps.resize(mix().nEnergyEqns());      
     }
 
     ~SetMixtureGSIO2gamma()
     {
         BOOST_TEST_MESSAGE("Destroying mixture");
         if(!mp_mix){
             delete mp_mix;
             mp_mix = NULL;
         }
     }
 
     Mutation::Mixture& mix() { return *mp_mix; }
 
     /**
      * Sets the mixture to an equilibrium state at the given T and P.  Assumes
      * that the variable set 1 of the state model accepts species densities and
      * temperature vectors.
      */
     void equilibrateMixture(double T, double P)
     {
 
         // Get equilibrium mole fractions
         mix().equilibriumComposition(T, P, &rhoi[0]);
 
         // Convert to species densities
         double rho = mix().density(T, P, &rhoi[0]);
         mix().convert<Mutation::Thermodynamics::X_TO_Y>(&rhoi[0], &rhoi[0]);
         for (int i = 0; i < mix().nSpecies(); ++i)
             rhoi[i] *= rho;
 
         tmps.assign(mix().nEnergyEqns(), T);
         mix().setState(&rhoi[0], &tmps[0], 1);
     }
 
     std::vector<double> getPartialDensities(){
         return rhoi;
     }
 
 private:
     Mutation::Mixture* mp_mix;
     std::vector<double> rhoi;
     std::vector<double> tmps;
 };

/*
 *
 */

class SetMixtureGSIair5gamma {
public:
    SetMixtureGSIair5gamma()
        : mp_mix(NULL)
    {
        BOOST_TEST_MESSAGE("Loading mixture for gsi test");
        mp_mix = new Mutation::Mixture("air5_gsi_gamma");

        rhoi.resize(mix().nSpecies());
        tmps.resize(mix().nEnergyEqns());      
    }

    ~SetMixtureGSIair5gamma()
    {
        BOOST_TEST_MESSAGE("Destroying mixture");
        if(!mp_mix){
            delete mp_mix;
            mp_mix = NULL;
        }
    }

    Mutation::Mixture& mix() { return *mp_mix; }

    /**
     * Sets the mixture to an equilibrium state at the given T and P.  Assumes
     * that the variable set 1 of the state model accepts species densities and
     * temperature vectors.
     */
    void equilibrateMixture(const double T, const double P)
    {
        // Get equilibrium mole fractions
        mix().equilibriumComposition(T, P, &rhoi[0]);

        // Convert to species densities
        double rho = mix().density(T, P, &rhoi[0]);
        mix().convert<Mutation::Thermodynamics::X_TO_Y>(&rhoi[0], &rhoi[0]);
        for (int i = 0; i < mix().nSpecies(); ++i)
            rhoi[i] *= rho;

        tmps.assign(mix().nEnergyEqns(), T);
        mix().setState(&rhoi[0], &tmps[0], 1);
    }

    std::vector<double> getPartialDensities(){
        return rhoi;
    }

private:
    Mutation::Mixture* mp_mix;
    std::vector<double> rhoi;
    std::vector<double> tmps;
};

//class SetMixtureGSI {
//public:
//    SetMixtureGSI()
//        : mp_mix(NULL)
//    {}
//
//    ~SetMixtureGSIair5gamma()
//    {
//        BOOST_TEST_MESSAGE("Destroying mixture");
//        if(!mp_mix)
//            delete mp_mix;
//    }
//
//    void initializeSetMixtureGSI(std::string mixture){
//        BOOST_TEST_MESSAGE("Loading mixture for gsi test");
//        mp_mix = new Mutation::Mixture(mixture);
//
//        rhoi.resize(mix().nSpecies());
//        tmps.resize(mix().nEnergyEqns());      
//    }
//
//    Mutation::Mixture& mix() { return *mp_mix; }
//
//    /**
//     * Sets the mixture to an equilibrium state at the given T and P.  Assumes
//     * that the variable set 1 of the state model accepts species densities and
//     * temperature vectors.
//     */
//    void equilibrateMixture(double T, double P)
//    {
//
//        // Get equilibrium mole fractions
//        mix().equilibriumComposition(T, P, &rhoi[0]);
//
//        // Convert to species densities
//        double rho = mix().density(T, P, &rhoi[0]);
//        mix().convert<Mutation::Thermodynamics::X_TO_Y>(&rhoi[0], &rhoi[0]);
//        for (int i = 0; i < mix().nSpecies(); ++i)
//            rhoi[i] *= rho;
//
//        tmps.assign(mix().nEnergyEqns(), T);
//        mix().setState(&rhoi[0], &tmps[0], 1);
//    }
//
//    std::vector<double> getPartialDensities(){
//        return rhoi;
//    }
//
//private:
//    Mutation::Mixture* mp_mix;
//    std::vector<double> rhoi;
//    std::vector<double> tmps;
//};
