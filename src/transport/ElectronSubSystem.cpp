/**
 * @file ElectronSubSystem.cpp
 *
 * @brief Implements ElectronSubSystem class.
 */

/*
 * Copyright 2016 von Karman Institute for Fluid Dynamics (VKI)
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

#include "ElectronSubSystem.h"

namespace Mutation {
    namespace Transport {

//==============================================================================

double  ElectronSubSystem::electricConductivity(int order)
{
    assert(order == 1 || order == 2);

    if (!m_thermo.hasElectrons())
        return 0.0;

    const double xe = m_thermo.X()[0];
    const double me = m_collisions.mass()(0);
    const double Te = m_thermo.Te();
    const double fac = 3.*xe*QE*QE/(16.*KB*Te)*std::sqrt(TWOPI*KB*Te/me);

    // First order
    if (order == 1)
        return fac / m_collisions.Lee<1>()(0,0);

    // Second order
    Eigen::Matrix2d L = m_collisions.Lee<2>();
    return fac / (L(0,0)-L(0,1)*L(0,1)/L(1,1));
}

//==============================================================================

Eigen::Vector3d ElectronSubSystem::electricConductivityB(int order)
{
    switch (order) {
    case 1: return sigmaB<1>();
    case 2: return sigmaB<2>();
    case 3: return sigmaB<3>();
    default:
        std::cout << "Warning: invalid order for sigma.  ";
        std::cout << "Using order 2..." << std::endl;
        return sigmaB<2>();
    }
}

//==============================================================================

double ElectronSubSystem::electronThermalConductivity(int order)
{
    assert(order == 2 || order == 3);

    if (!m_thermo.hasElectrons())
        return 0.0;

    const double xe = m_thermo.X()[0];
    const double Te = m_thermo.Te();
    const double fac = 75.*KB/64.*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));

    // 2nd order solution
    if (order == 2)
        return fac * xe / m_collisions.Lee<2>()(1,1);

    // 3rd order solution
    Eigen::Matrix3d L = m_collisions.Lee<3>();
    return fac * xe * L(2,2) / (L(1,1)*L(2,2) - L(1,2)*L(1,2));
}

//==============================================================================

template <>
Eigen::Vector3d ElectronSubSystem::electronThermalConductivityB<1>() {
    return Eigen::Vector3d::Zero();
}

//==============================================================================

Eigen::Vector3d ElectronSubSystem::electronThermalConductivityB(int order)
{
    switch (order) {
    case 1: return electronThermalConductivityB<1>();
    case 2: return electronThermalConductivityB<2>();
    case 3: return electronThermalConductivityB<3>();
    default:
        std::cout << "Warning: invalid order for electron thermal conductivity.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronThermalConductivityB<2>();
    }
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
