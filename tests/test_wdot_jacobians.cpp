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

//    SECTION("Air 5 mechanism jacobian wrt to density.") {
        // Initializing library
        Mixture mix("air5_jacobian_NASA-9_ChemNonEq1T");

        // Loop for multiple temperatures
        const double Tref = 300.;  // K

        const int ns = 5; //mix.nSpecies();
        const int nr = mix.nReactions();

        ArrayXd mw(ns); mw = mix.speciesMw();
        ArrayXd rhoi(ns); rhoi = mw; // One mole of everything.

        // ArrayXd jac_mpp(ns*ns);
        ArrayXd  jac_mpp(ns*ns);
        ArrayXXd jac(ns, ns);
        ArrayXXd dfrdrhoi(nr, ns);
        ArrayXXd dbrdrhoi(nr, ns);

        jac_mpp.setZero();

// Maybe loop. (Double loop to perturb everything, eg T, rhoi)
        mix.setState(rhoi.data(), &Tref, set_state_rhoi_T);
        mix.jacobianRho(jac_mpp.data());

        // Calculating jacobian wrt to rho for comparison
        ArrayXd kf(nr); ArrayXd kb(nr);
        mix.forwardRateCoefficients(kf.data());
        mix.backwardRateCoefficients(kb.data());

        dfrdrhoi.setZero(); dbrdrhoi.setZero();
        // Nitrogen dissociation: N2 + N = N + N + N
        // F: kf*[N2][N]
        dfrdrhoi(0,0) += kf(0)* rhoi(3) / (mw(3)*mw(0));   // wrt to N
        dfrdrhoi(0,3) += kf(0)* rhoi(0) / (mw(3)*mw(0));   // wrt to N2
        // B: kb*[N]*[N]*[N]
        dbrdrhoi(0,0) += kb(0)* 3*rhoi(0)*rhoi(0) / (mw(0)*mw(0)*mw(0));

        // Oxygen dissociation: O2 + M = O + O + M
        // F: kf*[O2][N]
        dfrdrhoi(1,0) += kf(1)* rhoi(4) / (mw(4)*mw(0)) * 5.;   // wrt to N
        dfrdrhoi(1,4) += kf(1)* rhoi(0) / (mw(4)*mw(0)) * 5.;   // wrt to O2
        // F: kf*[O2][O]
        dfrdrhoi(1,1) += kf(1)* rhoi(4) / (mw(4)*mw(1)) * 5.;   // wrt to O
        dfrdrhoi(1,4) += kf(1)* rhoi(1) / (mw(4)*mw(1)) * 5.;   // wrt to O2
        // F: kf*[O2][NO]
        dfrdrhoi(1,2) += kf(1)* rhoi(4) / (mw(4)*mw(2)) * 1.;   // wrt to NO
        dfrdrhoi(1,4) += kf(1)* rhoi(2) / (mw(4)*mw(2)) * 1.;   // wrt to O2
        // F: kf*[O2][N2]
        dfrdrhoi(1,3) += kf(1)* rhoi(4) / (mw(4)*mw(3)) * 1.;;   // wrt to N2
        dfrdrhoi(1,4) += kf(1)* rhoi(3) / (mw(4)*mw(3)) * 1.;;   // wrt to O2
        // F: kf*[O2][O2]
        dfrdrhoi(1,4) += kf(1)* 2*rhoi(2) / (mw(4)*mw(2)) * 1.;; // wrt to O2
        // B: kb*[O][O][N]
        dbrdrhoi(1,0) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(0)) * 5.;
        dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(0) / (mw(1)*mw(1)*mw(0)) * 5.;
        // B: kb*[O][O][O]
        dbrdrhoi(1,1) += kb(1)* 3*rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(1)) * 5.;
        // B: kb*[O][O][NO]
        dbrdrhoi(1,2) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(2)) * 1.;;
        dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(2) / (mw(1)*mw(1)*mw(2)) * 1.;;
        // B: kb*[O][O][N2]
        dbrdrhoi(1,3) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(3)) * 1.;;
        dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(3) / (mw(1)*mw(1)*mw(3)) * 1.;;
        // B: kb*[O][O][O2]
        dbrdrhoi(1,4) += kb(1)* rhoi(1)*rhoi(1) / (mw(1)*mw(1)*mw(4)) * 1.;;
        dbrdrhoi(1,1) += kb(1)* 2*rhoi(1)*rhoi(4) / (mw(1)*mw(1)*mw(4)) * 1.;;

        // Zeldovich: N2 + O = NO + N
        // F: kf[N2][O]
        dfrdrhoi(2,1) = kf(2)*rhoi(3) / (mw(3)*mw(1));
        dfrdrhoi(2,3) = kf(2)*rhoi(1) / (mw(3)*mw(1));
        // B: kb[NO][N]
        dbrdrhoi(2,0) = kb(2)*rhoi(2) / (mw(2)*mw(0));
        dbrdrhoi(2,2) = kb(2)*rhoi(0) / (mw(2)*mw(0));
        // Zeldovich: O + NO = O2 + N
        // F: kf[O][NO]
        // dfrdrhoi(3,1) = kf(3)*rhoi(2) / (mw(1)*mw(2));
        // dfrdrhoi(3,2) = kf(3)*rhoi(1) / (mw(1)*mw(2));
        // B: kb[O2][N]
        // dbrdrhoi(3,0) = kb(3)*rhoi(4) / (mw(4)*mw(0))
        // dbrdrhoi(3,4) = kb(3)*rhoi(0) / (mw(4)*mw(0))

        // Building the jacobian
        jac(0,0) = mw(0)*(
            +2*(dfrdrhoi(0,0) - dbrdrhoi(0,0))  // R1
            -1*(dbrdrhoi(2,0)));
        jac(0,1) = mw(0)*(
            -1*(dfrdrhoi(0,0)));
        jac(0,2) = mw(0)*(
            +1*(dbrdrhoi(2,2)));
        jac(0,3) = mw(0)*(
            +2*(dfrdrhoi(0,3))   // R1
            +1*(dfrdrhoi(2,3)));
        jac(0,4) = 0.;
        jac(1,0) = mw(1)*(
            +2*(dfrdrhoi(1,0) - dbrdrhoi(1,0)) // R2
            +1*(dbrdrhoi(2,0)));
        jac(1,1) = mw(1)*(
            +2*(dfrdrhoi(1,1) - dbrdrhoi(1,1))  // R2
            -1*(dfrdrhoi(2,1)));
        jac(1,2) = mw(1)*(
            +2*(dfrdrhoi(1,2) - dbrdrhoi(1,2)) // R2
            +1*(dbrdrhoi(2,2)));
        jac(1,3) = mw(1)*(
            +2*(dfrdrhoi(1,3) - dbrdrhoi(1,3))) // R2
            -1*(dfrdrhoi(2,3));
        jac(1,4) = mw(1)*(
            +2*(dfrdrhoi(1,4) - dbrdrhoi(1,4))); // R2
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
            -1*(dfrdrhoi(0,0) - dbrdrhoi(0,0))); // R1
        jac(3,1) = 0.;
        jac(3,2) = 0.;
        jac(3,3) = mw(3)*(
            -1*(dfrdrhoi(0,3))); // R1
        jac(3,4) = 0.;
        jac(4,0) = mw(4)*(
            -1*(dfrdrhoi(1,0) - dbrdrhoi(1,0))
            +1*(dbrdrhoi(2,0)));
        jac(4,1) = mw(4)*(
            -1*(dfrdrhoi(1,1) - dbrdrhoi(1,1))
            -1*(dfrdrhoi(2,1)));
        jac(4,2) = mw(4)*(
            -1*(dfrdrhoi(1,2) - dbrdrhoi(1,2))
            -1*(dfrdrhoi(2,3)));
        jac(4,3) = mw(4)*(
            -1*(dfrdrhoi(1,3) - dbrdrhoi(1,3))
            -1*(dfrdrhoi(2,3)));
        jac(4,4) = mw(4)*(
            -1*(dfrdrhoi(1,4) - dbrdrhoi(1,4)));

        std::cout << "Test" << std::endl;
        std::cout << jac << std::endl << std::endl;
        std::cout << "M++" << std::endl;
        for (int i = 0; i < ns; i++){
            for (int j = 0; j < ns; j++){
                std::cout << std::setw(12) << jac_mpp(i*ns+j) << " ";
            }
            std::cout << std::endl;
        }

        // Making sure everything is correct
        for (int i = 0; i < ns; i++)
            for (int j = 0; j < ns; j++)
                CHECK(jac(i,j) == Approx(jac_mpp(i*ns +j)).epsilon(tol));
}
