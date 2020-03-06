/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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

//==============================================================================

/// Helper class for diffusion matrix tests
class DiffusionMatrixTests
{
public:

    /**
     * Checks that species diffusion fluxes, computed by multiplying the
     * diffusion matrix times species densities, sum to zero.
     */
    static void fluxesSumToZero(Mixture& mix);
};

//==============================================================================

/**
 * Checks that species mass diffusion fluxes sum to zero when obtained by
 * multiplying the diffusion matrix by suitable driving forces.
 * \f[
 * \sum_i J_i = - \sum_i \rho_i \sum_j D_{ij} d'_j = 0
 * \f]
 * Driving forces are simulated as the mole fraction gradient in an equilibrium
 * mixture at constant pressure and elemental composition.  Under, these
 * conditions, the driving force is then
 * \f[
 * d'_j = \nabla x_j = \frac{\partial x_j}{\partial T} \nabla T,
 * \f]
 * where \f$\nabla T\f$ is chosen to make the maximum driving force equal to 1.
 */
TEST_CASE("DiffusionMatrix yields diffusion fluxes which sum to zero",
    "[transport][DiffusionMatrix]"
)
{
    MIXTURE_LOOP(
        //SECTION("Exact diffusion matrix") {
        //    mix.setDiffusionMatrixAlgo("Exact");
        //    fluxesSumToZero(mix);
        //}

        SECTION("Ramshaw diffusion matrix") {
            mix.setDiffusionMatrixAlgo("Ramshaw");
            DiffusionMatrixTests::fluxesSumToZero(mix);
        }
    )
}

//==============================================================================

void DiffusionMatrixTests::fluxesSumToZero(Mixture& mix)
{
    const int ns = mix.nSpecies();
    VectorXd dx(ns);
    VectorXd Ji(ns);

    const double tol = std::numeric_limits<double>::epsilon()*1000;

    // Test over several temperatures and pressures
    EQUILIBRATE_LOOP(
        // Compute temperature gradient as driving forces
        mix.dXidT(dx.data());

        // Scale to make the temperature gradient larger and make sure they
        // sum to zero
        dx *= 1.0 / dx.array().abs().maxCoeff();
        dx[0] -= dx.sum();

        // Compute the diffusion velocities
        Ji = -mix.diffusionMatrix() * dx; // <- diffusion velocities
        Ji.array() *= mix.density() * Map<const ArrayXd>(mix.Y(),ns);

        // Make sure the sum is zero
        INFO("Ji = " << Ji << "\ndx = " << dx);
        CHECK(Ji.sum() == Approx(0.0).margin(tol));
    )
}

//==============================================================================

