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

TEST_CASE("Surface chemical rates compared with analytical solution.", "[gsi]")
{
    const double tol = std::numeric_limits<double>::epsilon();
    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Catalysis Gamma Heteronuclear Reactions Single Temperature.")
    {
        // Setting up M++
        MixtureOptions opts("smb_air5_RRHO_ChemNonEq1T");
        Mixture mix(opts);
        const size_t ns = mix.nSpecies();
        VectorXd mol_mass = mix.speciesMw();

        const size_t nr = 3;
        CHECK(nr == mix.nSurfaceReactions());

        // Conditions Surface
        size_t nT = 6;
        for (int iT = 0; iT < nT; iT++) {
            double Tsurf = 300. + iT * 200;

            // Conditions dissociated
            double T = 3000.;
            double P = 100.;
            VectorXd rho(ns);
            mix.equilibrate(T, P);
            mix.densities(rho.data());

            // Catalysis
            const double gamma_N2 = .001;
            const double gamma_O2 = .001;
            const double gamma_NO = .002;
            const double gamma_ON = .003;

            double flux_N = sqrt((RU*Tsurf)/(2*PI*mol_mass(0)));
            double flux_O = sqrt((RU*Tsurf)/(2*PI*mol_mass(1)));
            VectorXd omega(ns);

            // Homonuclear reactions
            omega(0) = gamma_N2*flux_N*rho(0);
            omega(1) = gamma_O2*flux_O*rho(1);
            omega(2) = 0.;
            omega(3) = -gamma_N2*flux_N*rho(0);
            omega(4) = -gamma_O2*flux_O*rho(1);

            // Heteronuclear reactions
            if (gamma_NO * flux_N * rho(0) / mol_mass(0)
              < gamma_ON * flux_O * rho(1) / mol_mass(1)) {
                omega(0) += gamma_NO*flux_N*rho(0);
                omega(1) += gamma_NO*flux_N*rho(0)*mol_mass(1)/mol_mass(0);
                omega(2) -= gamma_NO*flux_N*rho(0)*mol_mass(2)/mol_mass(0);
            } else {
                omega(0) += gamma_ON*flux_O*rho(1)*mol_mass(0)/mol_mass(1);
                omega(1) += gamma_ON*flux_O*rho(1);
                omega(2) -= gamma_ON*flux_O*rho(1)*mol_mass(2)/mol_mass(1);
            }

            // Getting it from Mutation
            const size_t set_state_with_rhoi_T = 1;
            mix.setSurfaceState(rho.data(), &Tsurf, set_state_with_rhoi_T);
            VectorXd wdot(ns);
            mix.surfaceReactionRates(wdot.data());
            for (int i = 0; i < ns; i++) {
                INFO("The analytical solution for species " << i
                << " is equal to " << omega(i) << " while Mutation++ gives "
                << wdot(i) << ".\n");
                CHECK(wdot(i) == Approx(omega(i)).epsilon(tol));
            }
        }
    }

    SECTION("Catalysis Gamma Full Ion Recombination Reactions Single Temperature.")
    {
        // Setting up M++
        MixtureOptions opts("smb_air11_RRHO_ChemNonEq1T");
        Mixture mix(opts);
        const size_t ns = mix.nSpecies();
        VectorXd mol_mass = mix.speciesMw();

        const size_t nr = 8;
        CHECK(nr == mix.nSurfaceReactions());

        // Conditions Surface
         size_t nT = 6;
        for (int iT = 0; iT < nT; iT++) {
            double Tsurf = 300. + iT * 200;

            // Conditions dissociated
            double T = 3000.;
            double P = 100.;
            VectorXd rho(ns);
            mix.equilibrate(T, P);
            mix.densities(rho.data());

            // Catalysis
            const double gamma_N2 = .001;
            const double gamma_O2 = .001;
            const double gamma_NO = .002;
            const double gamma_ON = .003;
            const double gamma_em = 1.;

            double flux_N = sqrt((RU*Tsurf)/(2*PI*mol_mass(1)));
            double flux_O = sqrt((RU*Tsurf)/(2*PI*mol_mass(2)));
            VectorXd omega(ns);

            // Homonuclear reactions
            omega.setZero();
            omega(1) = gamma_N2*flux_N*rho(1);
            omega(2) = gamma_O2*flux_O*rho(2);
            omega(4) = -gamma_N2*flux_N*rho(1);
            omega(5) = -gamma_O2*flux_O*rho(2);

            // Heteronuclear reaction
            if (gamma_NO * flux_N * rho(1) / mol_mass(1)
              < gamma_ON * flux_O * rho(2) / mol_mass(2)) {
                omega(1) += gamma_NO*flux_N*rho(1);
                omega(2) += gamma_NO*flux_N*rho(1)*mol_mass(2)/mol_mass(1);
                omega(3) -= gamma_NO*flux_N*rho(1)*mol_mass(3)/mol_mass(1);
            } else {
                omega(1) += gamma_ON*flux_O*rho(2)*mol_mass(1)/mol_mass(2);
                omega(2) += gamma_ON*flux_O*rho(2);
                omega(3) -= gamma_ON*flux_O*rho(2)*mol_mass(3)/mol_mass(2);
            }

            double flux_em = sqrt((RU*Tsurf)/(2*PI*mol_mass(0)));
            double flux_Np = sqrt((RU*Tsurf)/(2*PI*mol_mass(6)));
            double flux_Op = sqrt((RU*Tsurf)/(2*PI*mol_mass(7)));
            double flux_NOp = sqrt((RU*Tsurf)/(2*PI*mol_mass(8)));
            double flux_N2p = sqrt((RU*Tsurf)/(2*PI*mol_mass(9)));
            double flux_O2p = sqrt((RU*Tsurf)/(2*PI*mol_mass(10)));

            CHECK(flux_Np * rho(6) / mol_mass(6) <
                  flux_em * rho(0) / mol_mass(0));
            CHECK(flux_Op * rho(7) / mol_mass(7) <
                  flux_em * rho(0) / mol_mass(0));
            CHECK(flux_NOp * rho(8) / mol_mass(8) <
                  flux_em * rho(0) / mol_mass(0));
            CHECK(flux_N2p * rho(9) / mol_mass(9) <
                  flux_em * rho(0) / mol_mass(0));
            CHECK(flux_O2p * rho(10) / mol_mass(10) <
                  flux_em * rho(0) / mol_mass(0));

            // Full Ion Recombination
            omega(1) -= gamma_em * flux_Np *
                        rho(6) * mol_mass(1)/mol_mass(6);
            omega(2) -= gamma_em * flux_Op *
                        rho(7) * mol_mass(2)/mol_mass(7);
            omega(3) -= gamma_em * flux_NOp *
                        rho(8) * mol_mass(3)/mol_mass(8);
            omega(4) -= gamma_em * flux_N2p *
                        rho(9) * mol_mass(4)/mol_mass(9);
            omega(5) -= gamma_em * flux_O2p *
                        rho(10) * mol_mass(5)/mol_mass(10);
            omega(6) += gamma_em * flux_Np * rho(6);
            omega(7) += gamma_em * flux_Op * rho(7);
            omega(8) += gamma_em * flux_NOp * rho(8);
            omega(9) += gamma_em * flux_N2p * rho(9);
            omega(10) += gamma_em * flux_O2p * rho(10);

            omega(0) += gamma_em*flux_Np*rho(6) * mol_mass(0)/mol_mass(6)
                      + gamma_em*flux_Op*rho(7) * mol_mass(0)/mol_mass(7)
                      + gamma_em*flux_NOp*rho(8) * mol_mass(0)/mol_mass(8)
                      + gamma_em*flux_N2p*rho(9) * mol_mass(0)/mol_mass(9)
                      + gamma_em*flux_O2p*rho(10) * mol_mass(0)/mol_mass(10);

            // Getting it from Mutation
            const size_t set_state_with_rhoi_T = 1;
            mix.setSurfaceState(rho.data(), &Tsurf, set_state_with_rhoi_T);
            VectorXd wdot(ns);
            mix.surfaceReactionRates(wdot.data());
            for (int i = 0; i < ns; i++) {
                INFO("The analytical solution for species " << i
                << " is equal to " << omega(i) << " while Mutation++ gives "
                << wdot(i) << ".\n");
                CHECK(wdot(i) == Approx(omega(i)).epsilon(tol));
            }
        }
    }

    SECTION("Phenomenological Ablation for Carbon Oxidation, Nitridation and Sublimation.")
    {
        // Setting up M++
        MixtureOptions opts("smb_aircarbon11_RRHO_ChemNonEq1T");
        Mixture mix(opts);
        const size_t ns = mix.nSpecies();
        VectorXd mol_mass = mix.speciesMw();

        const size_t nr = 4;
        CHECK(nr == mix.nSurfaceReactions());

        // Conditions Surface
        size_t nT = 1;
        for (int iT = 0; iT < nT; iT++) {
            double Tsurf = 500. + iT * 500;

            // Conditions dissociated
            double T = 3000.;
            double P = 100.;
            VectorXd rho(ns);
            mix.equilibrate(T, P);
            mix.densities(rho.data());

            // Catalysis
            const double gamma_CO = 0.63 * std::exp(-1160./Tsurf);
            const double gamma_CN = .003;
            const double gamma_N2 = .001;

            double flux_N = sqrt((RU*Tsurf)/(2*PI*mol_mass(0)));
            double flux_O = sqrt((RU*Tsurf)/(2*PI*mol_mass(1)));

            VectorXd omega(ns);
            omega.setZero();

            // Catalysis reactions
            omega(0) += gamma_N2*flux_N*rho(0);
            omega(3) -= gamma_N2*flux_N*rho(0);

            // Oxidation
            omega(1) += gamma_CO*flux_O*rho(1);
            omega(9) -= gamma_CO*flux_O*rho(1) * mol_mass(9)/mol_mass(1);

            // Nitridation
            omega(0) += gamma_CN*flux_N*rho(0);
            omega(8) -= gamma_CN*flux_N*rho(0) * mol_mass(8)/mol_mass(0);

            // Sublimation
            double vap_coef = 0.1;
            double Tact   = 90845.;
            double pre_exp = 5.19e15;
            double vtherm = sqrt(8*RU*Tsurf / (PI*mol_mass(7)));
            double p_sat = pre_exp * exp(-Tact/Tsurf);
            double rho_sat = p_sat * mol_mass(7) / (RU * Tsurf);
            omega(7) -= vap_coef * (rho_sat-rho(7)) * vtherm / 4.;

            // Getting it from Mutation
            const size_t set_state_with_rhoi_T = 1;
            mix.setSurfaceState(rho.data(), &Tsurf, set_state_with_rhoi_T);
            VectorXd wdot(ns);
            mix.surfaceReactionRates(wdot.data());
            for (int i = 0; i < ns; i++) {
                INFO("The analytical solution for species " << i
                << " is equal to " << omega(i) << " while Mutation++ gives "
                << wdot(i) << ".\n");
                CHECK(wdot(i) == Approx(omega(i)).epsilon(tol));
            }
        }
    }

    SECTION("Catalysis Gamma Homonuclear Reactions.")
    {
        // Setting up M++
        MixtureOptions opts("smb_o2_RRHO_ChemNonEq1T");
        Mixture mix(opts);
        const size_t ns = mix.nSpecies();
        double mol_mass = mix.speciesMw(0);

        const size_t nr = 1;
        CHECK(nr == mix.nSurfaceReactions());

        // Conditions Surface
        size_t nT = 4;
        for (int iT = 0; iT < nT; iT++) {
            double Tsurf = 300. + iT * 200;

            // Conditions dissociated
            double T = 3000.;
            double P = 100.;
            VectorXd rho(ns);
            mix.equilibrate(T, P);
            mix.densities(rho.data());

            // Catalysis
            const double gamma = 1.;
            double flux = sqrt((RU*Tsurf)/(2*PI*mol_mass));
            VectorXd omega(ns);
            omega(0) =  gamma*flux*rho(0);
            omega(1) = -gamma*flux*rho(0);

            // Getting it from Mutation
            const size_t set_state_with_rhoi_T = 1;
            mix.setSurfaceState(rho.data(),&Tsurf, set_state_with_rhoi_T);
            VectorXd wdot(ns);
            mix.surfaceReactionRates(wdot.data());
            for (int i = 0; i < ns; i++){
                INFO("The analytical solution for species " << i
                << " is equal to " << omega(i) << " while Mutation++ gives "
                << wdot(i) << ".\n");
                CHECK(wdot(i) == Approx(omega(i)).epsilon(tol));
            }
        }
    }

}
