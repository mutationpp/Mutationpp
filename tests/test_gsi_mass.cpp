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

TEST_CASE("Solution of the MassBalanceSolver is converged.",
    "[gsi]"
)
{
    const double tol = 100*std::numeric_limits<double>::epsilon();

    MIXTURE_GSI_MASS_LOOP
    (
        // Setting up
        const size_t set_state_with_rhoi_T = 1;
        const size_t pos_T_trans = 0;
        size_t ns = mix.nSpecies();
        size_t nT = mix.nEnergyEqns();

        // Conditions T = 3000K and p = 100Pa
        VectorXd Teq = VectorXd::Constant(nT, 3000.);
        double Peq = 100.; // Pa
        mix.equilibrate(Teq(pos_T_trans), Peq);
        VectorXd rhoi_s(ns);
        mix.densities(rhoi_s.data());

        VectorXd T_s = VectorXd::Constant(nT, 300.);

        // Setting number of iterations
        const int iter = 100;
        mix.setIterationsSurfaceBalance(iter);

        // Mass gradient
        VectorXd xi_e(ns);
        xi_e = Map<const VectorXd>(mix.X(), ns);
        double dx = 1.;

        // Solving mass balance
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        mix.setDiffusionModel(xi_e.data(), dx);

        mix.solveSurfaceBalance();
        mix.getSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

        // Verifying the solution gives low residual in the balance equations
        mix.setState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        VectorXd xi_s(ns);
        xi_s = Map<const VectorXd>(mix.X(), ns);

        // Compute diffusion velocities
        VectorXd dxidx(ns);
        dxidx = (xi_s - xi_e) / dx;
        VectorXd vdi(ns);
        double E = 0.;
        mix.stefanMaxwell(dxidx.data(), vdi.data(), E);

        // Get surface production rates
        VectorXd wdot(ns);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        mix.surfaceReactionRates(wdot.data());

        // Blowing flux (should be zero for catalysis)
        double mblow;
        mix.getMassBlowingRate(mblow);

        // Building mass balance functions
        VectorXd F(ns);
        double rho = rhoi_s.sum();
        F = (rhoi_s/rhoi_s.sum())*mblow + rhoi_s.cwiseProduct(vdi) - wdot;

        // Compute error
         double err = F.lpNorm<Infinity>();

        CHECK(err == Approx(0.0).epsilon(tol));
    )
}
