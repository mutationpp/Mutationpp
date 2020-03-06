/**
 * @file ElectronSubSystem.cpp
 *
 * @brief Implements ElectronSubSystem class.
 */

/*
 * Copyright 2016-2020 von Karman Institute for Fluid Dynamics (VKI)
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

double ElectronSubSystem::electricConductivity(int order)
{
    if (!m_thermo.hasElectrons())
        return 0.0;

    const double fac =
        m_thermo.numberDensity()*m_thermo.X()[0]*QE*QE/(KB*m_thermo.Te());

    return fac * electronDiffusionCoefficient(order);
}

//==============================================================================

Eigen::Vector3d ElectronSubSystem::electricConductivityB(int order)
{
    if (!m_thermo.hasElectrons())
        return Eigen::Vector3d::Zero();

    const double fac =
        m_thermo.numberDensity()*m_thermo.X()[0]*QE*QE/(KB*m_thermo.Te());

    return fac * electronDiffusionCoefficientB(order);
}

//==============================================================================

double ElectronSubSystem::electronDiffusionCoefficient(int order)
{
    switch (order) {
    case 1: return electronDiffusionCoefficient<1>();
    case 2: return electronDiffusionCoefficient<2>();
    case 3: return electronDiffusionCoefficient<3>();
    default:
        std::cout << "Warning: invalid order for electron diffusion coefficient.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronDiffusionCoefficient<3>();
    }
}

//==============================================================================

Eigen::Vector3d ElectronSubSystem::electronDiffusionCoefficientB(int order)
{
    switch (order) {
    case 1: return electronDiffusionCoefficientB<1>();
    case 2: return electronDiffusionCoefficientB<2>();
    case 3: return electronDiffusionCoefficientB<3>();
    default:
        std::cout << "Warning: invalid order for electron diffusion coefficient.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronDiffusionCoefficientB<3>();
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
        return electronThermalConductivityB<3>();
    }
}

//==============================================================================

const Eigen::VectorXd& ElectronSubSystem::alpha(int order)
{
    switch (order) {
    case 1: return alpha<1>();
    case 2: return alpha<2>();
    case 3: return alpha<3>();
    default:
        std::cout << "Warning: invalid order for alpha coefficients.  ";
        std::cout << "Using order 3..." << std::endl;
        return alpha<3>();
    }
}

//==============================================================================

const Eigen::Matrix<double,-1,3>& ElectronSubSystem::alphaB(int order)
{
    switch (order) {
    case 1: return alphaB<1>();
    case 2: return alphaB<2>();
    case 3: return alphaB<3>();
    default:
        std::cout << "Warning: invalid order for alpha coefficients.  ";
        std::cout << "Using order 3..." << std::endl;
        return alphaB<3>();
    }
}

//==============================================================================

template <>
double ElectronSubSystem::electronThermalDiffusionRatio<1>() {
    return 0.0;
}

//==============================================================================

double ElectronSubSystem::electronThermalDiffusionRatio(int order)
{
    switch (order) {
    case 1: return electronThermalDiffusionRatio<1>();
    case 2: return electronThermalDiffusionRatio<2>();
    case 3: return electronThermalDiffusionRatio<3>();
    default:
        std::cout << "Warning: invalid order for electron thermal diffusion ratio.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronThermalDiffusionRatio<3>();
    }
}

//==============================================================================

template <>
Eigen::Vector3d ElectronSubSystem::electronThermalDiffusionRatioB<1>() {
    return Eigen::Vector3d::Zero();
}

//==============================================================================

Eigen::Vector3d ElectronSubSystem::electronThermalDiffusionRatioB(int order)
{
    switch (order) {
    case 1: return electronThermalDiffusionRatioB<1>();
    case 2: return electronThermalDiffusionRatioB<2>();
    case 3: return electronThermalDiffusionRatioB<3>();
    default:
        std::cout << "Warning: invalid order for electron thermal diffusion ratio.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronThermalDiffusionRatioB<3>();
    }
}

//==============================================================================

const Eigen::VectorXd& ElectronSubSystem::electronThermalDiffusionRatios2(int order)
{
    switch (order) {
    case 1: return electronThermalDiffusionRatios2<1>();
    case 2: return electronThermalDiffusionRatios2<2>();
    case 3: return electronThermalDiffusionRatios2<3>();
    default:
        std::cout << "Warning: invalid order for 2nd order electron thermal diffusion ratios.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronThermalDiffusionRatios2<3>();
    }
}

//==============================================================================

const Eigen::Matrix<double,-1,3>& ElectronSubSystem::electronThermalDiffusionRatios2B(int order)
{
    switch (order) {
    case 1: return electronThermalDiffusionRatios2B<1>();
    case 2: return electronThermalDiffusionRatios2B<2>();
    case 3: return electronThermalDiffusionRatios2B<3>();
    default:
        std::cout << "Warning: invalid order for 2nd order electron thermal diffusion ratios.  ";
        std::cout << "Using order 3..." << std::endl;
        return electronThermalDiffusionRatios2B<3>();
    }
}

//==============================================================================

double ElectronSubSystem::electronDiffusionCoefficient2(
    const Eigen::Ref<const Eigen::MatrixXd>& Dij, int order)
{
    const Eigen::VectorXd a = (*this).alpha(order);
    return (Dij * a).dot(a);
}

//==============================================================================

Eigen::Vector3d ElectronSubSystem::electronDiffusionCoefficient2B(
    const Eigen::Ref<const Eigen::MatrixXd>& Dij, int order)
{
    const Eigen::Matrix<double,-1,3> a = (*this).alphaB(order);
    Eigen::Vector3d Dee;

    Dee(0) = (Dij * a.col(0)).dot(a.col(0));

    Dee(1) = (Dij * a.col(1)).dot(a.col(1)) -
             (Dij * a.col(2)).dot(a.col(2));

    Dee(2) = (Dij * a.col(1)).dot(a.col(2)) +
             (Dij * a.col(2)).dot(a.col(1));

    return Dee;
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
