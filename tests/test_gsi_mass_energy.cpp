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

TEST_CASE("Solution of the MassEnergyBalanceSolver is converged.", "[gsi]")
{
    const double tol = std::numeric_limits<double>::epsilon();
    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Inert Surface in Radiative Equilibrium.")
    {
        // Mixture
        MixtureOptions opts("seb_inert_NASA9_ChemNonEq1T");
        Mixture mix(opts);

        // Setting up
        const size_t set_state_with_rhoi_T = 1;
        const size_t pos_T_trans = 0;
        size_t ns = mix.nSpecies();
        size_t nT = mix.nEnergyEqns();
        size_t neq = ns + nT;

        // Conditions T = 3000K and p = 100Pa
        VectorXd Teq = VectorXd::Constant(nT, 3000.);
        double Peq = 100.; // Pa
        mix.equilibrate(Teq(pos_T_trans), Peq);

        // Setting number of iterations for the solver
        const int iter = 100;
        mix.setIterationsSurfaceBalance(iter);

        // Mass gradient
        VectorXd xi_e(ns);
        xi_e = Map<const VectorXd>(mix.X(), ns);
        double dx = 1.e-3;
        mix.setDiffusionModel(xi_e.data(), dx);

        // Temperature gradient
        VectorXd T_e = Teq;
        mix.setGasFourierHeatFluxModel(T_e.data(), dx);

        // Initial conditions of the surface are the ones in the first
        // physical cell
        VectorXd rhoi_s(ns);
        mix.densities(rhoi_s.data());
        VectorXd T_s = VectorXd::Constant(nT, 1800.);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

        // Solve balance and request solution
        mix.solveSurfaceBalance();
        mix.getSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

        // Verifying the solution gives low residual in the balance equations
        mix.setState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        VectorXd xi_s(ns);
        xi_s = Map<const VectorXd>(mix.X(), ns);

        // Compute diffusion velocities
        VectorXd dxidx(ns);
        dxidx = (xi_s - xi_e) / dx;
        VectorXd vdi = VectorXd::Zero(ns);
        double E = 0.;
        mix.stefanMaxwell(dxidx.data(), vdi.data(), E);

        // Conductive heat flux
        VectorXd dTdx(nT);
        dTdx = (T_s - T_e) / dx;
        VectorXd lambda(nT);
        mix.frozenThermalConductivityVector(lambda.data());

        // Get surface production rates
        VectorXd wdot(ns);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        mix.surfaceReactionRates(wdot.data());

        // Surface radiation
        const double sigma = 2.*pow(PI, 5)*pow(KB, 4)/(15*pow(C0, 2)*pow(HP, 3));
        const double eps = 1.;
        VectorXd q_srad = VectorXd::Zero(nT) ;
        q_srad(pos_T_trans) = sigma * eps * pow(T_s(pos_T_trans), 4);

        // Building balance functions
        VectorXd F(neq);
        double rho = rhoi_s.sum();
        F.head(ns) = rhoi_s.cwiseProduct(vdi) - wdot;
        F.tail(nT) = lambda.cwiseProduct(dTdx) + q_srad;

        // Compute error
        double err = F.lpNorm<Infinity>();
        CHECK(err == Approx(0.0).margin(tol));
    }

    SECTION("Mass and Energy Balance Catalysis in Radiative Equilibrium.")
    {
        // Mixture
        MixtureOptions opts("seb_catalysis_NASA9_ChemNonEq1T");
        Mixture mix(opts);

        // Setting up
        const size_t set_state_with_rhoi_T = 1;
        const size_t pos_T_trans = 0;
        size_t ns = mix.nSpecies();
        size_t nT = mix.nEnergyEqns();
        size_t neq = ns + nT;

        // Conditions T = 3000K and p = 100Pa
        VectorXd Teq = VectorXd::Constant(nT, 3000.);
        double Peq = 100.; // Pa
        mix.equilibrate(Teq(pos_T_trans), Peq);

        // Setting number of iterations for the solver
        const int iter = 100;
        mix.setIterationsSurfaceBalance(iter);

        // Mass gradient
        VectorXd xi_e(ns);
        xi_e = Map<const VectorXd>(mix.X(), ns);
        double dx = 1.e-3;
        mix.setDiffusionModel(xi_e.data(), dx);

        // Temperature gradient
        VectorXd T_e = Teq;
        mix.setGasFourierHeatFluxModel(T_e.data(), dx);

        // Initial conditions of the surface are the ones in the first
        // physical cell
        VectorXd rhoi_s(ns);
        mix.densities(rhoi_s.data());
        VectorXd T_s = VectorXd::Constant(nT, 1800.);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

        // Solve balance and request solution
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

        // Conductive heat flux
        VectorXd dTdx(nT);
        dTdx = (T_s - T_e) / dx;
        VectorXd lambda(nT);
        mix.frozenThermalConductivityVector(lambda.data());

        // Get surface production rates
        VectorXd wdot(ns);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        mix.surfaceReactionRates(wdot.data());

        // Species and mixture enthalpies
        VectorXd v_hi(ns*nT);
        mix.getEnthalpiesMass(v_hi.data());
        VectorXd v_h(nT);
        v_h(pos_T_trans) = rhoi_s.dot(v_hi.head(ns));

        // Surface radiation
        const double sigma = 2.*pow(PI, 5)*pow(KB, 4)/(15*pow(C0, 2)*pow(HP, 3));
        const double eps = .86;
        VectorXd q_srad = VectorXd::Zero(nT);
        q_srad(pos_T_trans) = sigma*eps*pow(T_s(pos_T_trans), 4);

        // Chemical Energy Contribution
        VectorXd v_hi_rhoi_vi= VectorXd::Zero(nT);
        v_hi_rhoi_vi(pos_T_trans) = -v_hi.head(ns).dot(rhoi_s.cwiseProduct(vdi));

        // Building balance functions
        VectorXd F(neq);
        double rho = rhoi_s.sum();
        F.head(ns) = rhoi_s.cwiseProduct(vdi) - wdot;
        F.tail(nT) = lambda.cwiseProduct(dTdx) + q_srad + v_hi_rhoi_vi;

        // Compute error
        double err = F.lpNorm<Infinity>();
        CHECK(err == Approx(0.0).margin(tol));
    }

    SECTION("Mass and Energy Balance Ablation.")
    {
        // Mixture
        MixtureOptions opts("seb_oxidation_NASA9_ChemNonEq1T");
        Mixture mix(opts);

        // Setting up
        const size_t set_state_with_rhoi_T = 1;
        const size_t pos_T_trans = 0;
        size_t ns = mix.nSpecies();
        size_t nT = mix.nEnergyEqns();
        size_t neq = ns + nT;

        // Conditions T = 3000K and p = 100Pa
        VectorXd Teq = VectorXd::Constant(nT, 3000.);
        double Peq = 100.; // Pa
        mix.equilibrate(Teq(pos_T_trans), Peq);

        // Setting number of iterations for the solver
        const int iter = 100;
        mix.setIterationsSurfaceBalance(iter);

        // Mass gradient
        VectorXd xi_e(ns);
        xi_e = Map<const VectorXd>(mix.X(), ns);
        double dx = 1.e-3;
        mix.setDiffusionModel(xi_e.data(), dx);

        // Temperature gradient
        VectorXd T_e = Teq;
        mix.setGasFourierHeatFluxModel(T_e.data(), dx);

        // Initial conditions of the surface are the ones in the first
        // physical cell
        VectorXd rhoi_s(ns);
        mix.densities(rhoi_s.data());
        VectorXd T_s = VectorXd::Constant(nT, 1800.);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

        // Solve balance and request solution
        mix.solveSurfaceBalance();
        mix.getSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        double rho = rhoi_s.sum();

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

        // Conductive heat flux
        VectorXd dTdx(nT);
        dTdx = (T_s - T_e) / dx;
        VectorXd lambda(nT);
        mix.frozenThermalConductivityVector(lambda.data());

        // Get surface production rates
        VectorXd wdot(ns);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        mix.surfaceReactionRates(wdot.data());

        // Blowing flux (should be zero for catalysis)
        double mblow;
        mix.getMassBlowingRate(mblow);

        // Species and mixture enthalpies
        VectorXd v_hi(ns*nT);
        mix.getEnthalpiesMass(v_hi.data());
        VectorXd v_h(nT);
        v_h(pos_T_trans) = (rhoi_s/rho).dot(v_hi.head(ns));

        // Surface radiation
        const double sigma = 2.*pow(PI, 5)*pow(KB, 4)/(15*pow(C0, 2)*pow(HP, 3));
        const double eps = .86;
        VectorXd q_srad = VectorXd::Zero(nT);
        q_srad(pos_T_trans) = sigma*eps*pow(T_s(pos_T_trans), 4);

        // Chemical Energy Contribution
        VectorXd v_hi_rhoi_vi= VectorXd::Zero(nT);
        v_hi_rhoi_vi(pos_T_trans) = -v_hi.head(ns).dot(rhoi_s.cwiseProduct(vdi));

        // Building balance functions
        VectorXd F(neq);
        F.head(ns) = (rhoi_s/rho)*mblow + rhoi_s.cwiseProduct(vdi) - wdot;
        F.tail(nT) = -lambda.cwiseProduct(dTdx) - q_srad + mblow*v_h -
            v_hi_rhoi_vi;

        // Compute error
        double err = F.lpNorm<Infinity>();
        CHECK(err == Approx(0.0).margin(tol));
    }

    SECTION("Mass and Energy Balance with Two Temperatures at Equilibrium."){
        // Mixture
        MixtureOptions opts("seb_aircarbon11_ablation_NASA9_ChemNonEqTTv");
        Mixture mix(opts);

        // Setting up
        const size_t set_state_with_rhoi_T = 1;
        const size_t pos_T_trans = 0;
        const size_t pos_T_vib = 1;
        size_t ns = mix.nSpecies();
        size_t nT = mix.nEnergyEqns();
        size_t nE = 1;
        size_t neq = ns + nE;

        // Conditions T = 3000K and p = 100Pa
        VectorXd Teq = VectorXd::Constant(nT, 3000.);
        double Peq = 100.; // Pa
        mix.equilibrate(Teq(pos_T_trans), Peq);

        // Setting number of iterations for the solver
        const int iter = 100;
        mix.setIterationsSurfaceBalance(iter);

        // Mass gradient
        VectorXd xi_e(ns);
        xi_e = Map<const VectorXd>(mix.X(), ns);
        double dx = 1.e-3;
        mix.setDiffusionModel(xi_e.data(), dx);

        // Temperature gradient
        VectorXd T_e = Teq;
        mix.setGasFourierHeatFluxModel(T_e.data(), dx);

        // Initial conditions of the surface are the ones in the first
        // physical cell
        VectorXd rhoi_s(ns);
        mix.densities(rhoi_s.data());
        VectorXd T_s(nT);
        T_s(pos_T_trans) = 1800.;
        T_s(pos_T_vib)   = 2200.;
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

        // Solve balance and request solution
        mix.solveSurfaceBalance();
        mix.getSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        double rho = rhoi_s.sum();

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

        // Conductive heat flux
        VectorXd dTdx(nT);
        dTdx = (T_s - T_e) / dx;
        VectorXd lambda(nT);
        mix.frozenThermalConductivityVector(lambda.data());

        // Get surface production rates
        VectorXd wdot(ns);
        mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
        mix.surfaceReactionRates(wdot.data());

        // Blowing flux (should be zero for catalysis)
        double mblow;
        mix.getMassBlowingRate(mblow);

        // Species and mixture enthalpies
        VectorXd v_hi(ns*nT);
        mix.getEnthalpiesMass(v_hi.data());
        double h = (rhoi_s/rho).dot(v_hi.head(ns));

        // Surface radiation
        const double sigma = 2.*pow(PI, 5)*pow(KB, 4)/(15*pow(C0, 2)*pow(HP, 3));
        const double eps = .86;
        double q_srad = sigma*eps*pow(T_s(pos_T_trans), 4);

        // Chemical Energy Contribution
        double v_hi_rhoi_vi = -v_hi.head(ns).dot(rhoi_s.cwiseProduct(vdi));

        // Building balance functions
        VectorXd F(neq);
        F.head(ns) = (rhoi_s/rho)*mblow + rhoi_s.cwiseProduct(vdi) - wdot;
        F(nE) = -lambda.dot(dTdx) - q_srad + mblow*h - v_hi_rhoi_vi;
        F.tail(nT-1).setConstant(0.);

        // Compute error
        double err = F.lpNorm<Infinity>();
        CHECK(err == Approx(0.0).margin(tol));
    }
}
