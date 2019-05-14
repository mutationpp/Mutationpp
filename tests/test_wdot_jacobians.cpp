/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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

TEST_CASE("Chemical rates jacobians.", "[kinetics]")
{
    const double tol = 100*std::numeric_limits<double>::epsilon();
    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    const int set_state_rhoi_T = 1.;

    /* SECTION("Air 5 Jacobian rhoi") {
        Mixture mix("air5_jacobian_NASA-9_ChemNonEq1T");

        const int ns = mix.nSpecies();
        const int nr = mix.nReactions();

        ArrayXd mw(ns); mw = mix.speciesMw();
        ArrayXd rhoi(ns); rhoi = mw;

        ArrayXd  jac_mpp(ns*ns);
        ArrayXXd jac(ns, ns);
        ArrayXXd dfrdrhoi(3, ns);
        ArrayXXd dbrdrhoi(3, ns);

        jac_mpp.setZero();

        // Loop for multiple temperatures
        const int nT = 4;
        for (int iT = 0; iT < nT; iT++) {
            double T = 500*iT + 300.;  // K

            mix.setState(rhoi.data(), &T, set_state_rhoi_T);
            mix.jacobianRho(jac_mpp.data());

            // Calculating jacobian wrt to rho for comparison
            ArrayXd kf(nr); ArrayXd kb(nr);
            mix.forwardRateCoefficients(kf.data());
            mix.backwardRateCoefficients(kb.data());

            dfrdrhoi.setZero(); dbrdrhoi.setZero();

            // Nitrogen dissociation: N2 + N = N + N + N
            // F: kf*[N2][N]
            dfrdrhoi(0,0) += kf(0)* rhoi(3) / (mw(3)*mw(0));
            dfrdrhoi(0,3) += kf(0)* rhoi(0) / (mw(3)*mw(0));
            // B: kb*[N]*[N]*[N]
            dbrdrhoi(0,0) += kb(0)* 3*rhoi(0)*rhoi(0) / (mw(0)*mw(0)*mw(0));

            // Oxygen dissociation: O2 + M = O + O + M
            // F: kf*[O2][N]
            dfrdrhoi(1,0) += kf(1)* rhoi(4) / (mw(4)*mw(0)) * 5.;
            dfrdrhoi(1,4) += kf(1)* rhoi(0) / (mw(4)*mw(0)) * 5.;
            // F: kf*[O2][O]
            dfrdrhoi(1,1) += kf(1)* rhoi(4) / (mw(4)*mw(1)) * 5.;
            dfrdrhoi(1,4) += kf(1)* rhoi(1) / (mw(4)*mw(1)) * 5.;
            // F: kf*[O2][NO]
            dfrdrhoi(1,2) += kf(1)* rhoi(4) / (mw(4)*mw(2)) * 1.;
            dfrdrhoi(1,4) += kf(1)* rhoi(2) / (mw(4)*mw(2)) * 1.;
            // F: kf*[O2][N2]
            dfrdrhoi(1,3) += kf(1)* rhoi(4) / (mw(4)*mw(3)) * 1.;
            dfrdrhoi(1,4) += kf(1)* rhoi(3) / (mw(4)*mw(3)) * 1.;
            // F: kf*[O2][O2]
            dfrdrhoi(1,4) += kf(1)* 2*rhoi(2) / (mw(4)*mw(2)) * 1.;
            // B: kb*[O][O][N]
            dbrdrhoi(1,0) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(0)) * 5.;
            dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(0) / (mw(1)*mw(1)*mw(0)) * 5.;
            // B: kb*[O][O][O]
            dbrdrhoi(1,1) += kb(1)* 3*rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(1)) * 5.;
            // B: kb*[O][O][NO]
            dbrdrhoi(1,2) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(2)) * 1.;
            dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(2) / (mw(1)*mw(1)*mw(2)) * 1.;
            // B: kb*[O][O][N2]
            dbrdrhoi(1,3) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(3)) * 1.;
            dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(3) / (mw(1)*mw(1)*mw(3)) * 1.;
            // B: kb*[O][O][O2]
            dbrdrhoi(1,4) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(4)) * 1.;
            dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(4) / (mw(1)*mw(1)*mw(4)) * 1.;

            // Zeldovich: N2 + O = NO + N
            // F: kf[N2][O]
            dfrdrhoi(2,1) += kf(2)*rhoi(3) / (mw(3)*mw(1));
            dfrdrhoi(2,3) += kf(2)*rhoi(1) / (mw(3)*mw(1));
            // B: kb[NO][N]
            dbrdrhoi(2,0) += kb(2)*rhoi(2) / (mw(2)*mw(0));
            dbrdrhoi(2,2) += kb(2)*rhoi(0) / (mw(2)*mw(0));

            // Building the jacobian
            jac(0,0) = mw(0)*(
                +2*(dfrdrhoi(0,0) - dbrdrhoi(0,0))
                -1*(dbrdrhoi(2,0)));
            jac(0,1) = mw(0)*(
                +1*(dfrdrhoi(2,1)));
            jac(0,2) = mw(0)*(
                -1*(dbrdrhoi(2,2)));
            jac(0,3) = mw(0)*(
                +2*(dfrdrhoi(0,3))
                +1*(dfrdrhoi(2,3)));
            jac(0,4) = 0.;
            jac(1,0) = mw(1)*(
                +2*(dfrdrhoi(1,0) - dbrdrhoi(1,0))
                +1*(dbrdrhoi(2,0)));
            jac(1,1) = mw(1)*(
                +2*(dfrdrhoi(1,1) - dbrdrhoi(1,1))
                -1*(dfrdrhoi(2,1)));
            jac(1,2) = mw(1)*(
                +2*(dfrdrhoi(1,2) - dbrdrhoi(1,2))
                +1*(dbrdrhoi(2,2)));
            jac(1,3) = mw(1)*(
                +2*(dfrdrhoi(1,3) - dbrdrhoi(1,3))
                -1*(dfrdrhoi(2,3)));
            jac(1,4) = mw(1)*(
                +2*(dfrdrhoi(1,4) - dbrdrhoi(1,4)));
            jac(2,0) = mw(2)*(
                -1*(dbrdrhoi(2,0)));
            jac(2,1) = mw(2)*(
                +1*(dfrdrhoi(2,1)));
            jac(2,2) = mw(2)*(
                -1*(dbrdrhoi(2,2)));
            jac(2,3) = mw(2)*(
                +1*(dfrdrhoi(2,3)));
            jac(2,4) = 0;
            jac(3,0) = mw(3)*(
                -1*(dfrdrhoi(0,0) - dbrdrhoi(0,0))
                +1*(dbrdrhoi(2,0)));
            jac(3,1) = mw(3)*(
                -1*(dfrdrhoi(2,1)));
            jac(3,2) = mw(3)*(
                +1*(dbrdrhoi(2,2)));
            jac(3,3) = mw(3)*(
                -1*(dfrdrhoi(0,3)) // R1
                -1*(dfrdrhoi(2,3)));
            jac(3,4) = 0.;
            jac(4,0) = mw(4)*(
                -1*(dfrdrhoi(1,0) - dbrdrhoi(1,0))
                +0);
            jac(4,1) = mw(4)*(
                -1*(dfrdrhoi(1,1) - dbrdrhoi(1,1))
                -0);
            jac(4,2) = mw(4)*(
                -1*(dfrdrhoi(1,2) - dbrdrhoi(1,2))
                -0);
            jac(4,3) = mw(4)*(
                -1*(dfrdrhoi(1,3) - dbrdrhoi(1,3))
                -0);
            jac(4,4) = mw(4)*(
                -1*(dfrdrhoi(1,4)-dbrdrhoi(1,4)));

            // Making sure everything is correct
            for (int i = 0; i < ns; i++)
                for (int j = 0; j < ns; j++)
                    CHECK(jac(i,j) == Approx(jac_mpp(i*ns +j)).epsilon(tol));
        }
    } */

    SECTION("Air 5 Jacobian T") {
        Mixture mix("air5_test");

        const int ns = mix.nSpecies();
        const int nr = mix.nReactions();

        // Reaction coefficients
        double A = 7.0E+15;
        double n = -1.6;
        double Ta = 113200.;

        double T = 3000.;
        double invT = 1./T;

        ArrayXd mw(ns); mw = mix.speciesMw();
        ArrayXd rhoi(ns); rhoi = mw;
        ArrayXd gi(ns); ArrayXd hi(ns);

        mix.setState(rhoi.data(), &T, set_state_rhoi_T);

        ArrayXd vkf(nr); ArrayXd vkb(nr); ArrayXd vKeq(nr);
        ArrayXd vdkfdT(nr); ArrayXd vdkbdT(nr);

        mix.speciesSTGOverRT(T, gi.data());
        mix.speciesHOverRT(hi.data());

        vkf(0) = A * pow(T, n) * exp(-Ta/T);
        vdkfdT(0) = vkf(0)*invT*(n + invT*Ta);
        vKeq(0) = ONEATM/(RU*T)*exp(-(2*gi(0)-gi(3)));
        vkb(0) = vkf(0)/vKeq(0);
        double dlnkeqdT = 2*(hi(0) - 1) - (hi(3) - 1);
        vdkbdT(0) = 1/vKeq(0)*(vdkfdT(0) - vkf(0)/T * dlnkeqdT);

        ArrayXd kfmpp(nr); ArrayXd kbmpp(nr);
        ArrayXd dkfdTmpp(nr); ArrayXd dkbdTmpp(nr);

        mix.forwardRateCoefficients(kfmpp.data());
        mix.dkfdT(dkfdTmpp.data());
        mix.backwardRateCoefficients(kbmpp.data());
        mix.dkbdT(dkbdTmpp.data());

        // std::cout << "Kb = "    << vkb    << " " << kbmpp << std::endl;
        // std::cout << "Kf = "    << vkf    << " " << kfmpp << std::endl;
        // std::cout << "dkfdT = " << vdkfdT << " " << dkfdTmpp << std::endl;
        // std::cout << "dkbdT = " << vdkbdT << " " << dkbdTmpp << std::endl;

        CHECK(vkf(0) == Approx(kfmpp(0)).epsilon(tol));
        CHECK(vkb(0) == Approx(kbmpp(0)).epsilon(tol));
        CHECK(vdkfdT(0) == Approx(dkfdTmpp(0)).epsilon(tol));
        CHECK(vdkbdT(0) == Approx(dkbdTmpp(0)).epsilon(tol));
    }
}
