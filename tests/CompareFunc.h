/**
 * @file CompareFunc.h
 *
 * @brief Provides self registering function comparison classes.
 */

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

#ifndef TESTS_COMPARE_FUNC_H
#define TESTS_COMPARE_FUNC_H

#include "mutation++.h"
#include <Eigen/Dense>
#include <iostream>

/**
 * Base class for all comparison functions.
 */
class CompareFunc
{
public:
    typedef double ARGS;

    /// Returns name of this type.
    static std::string typeName() { return "CompareFunc"; }

    CompareFunc(double tol) : m_tolerance(tol) {}

    virtual ~CompareFunc() {}

    bool compare(Mutation::Mixture& mix, Eigen::VectorXd result)
    {
        Eigen::VectorXd v(result.size());
        compute(mix, v);
        bool matches =  v.isApprox(result, m_tolerance);

        if (!matches) {
            std::cout << "Comparison between {"
                 << v.transpose() << "} and {"
                 << result.transpose() << "} failed."
                 << std::endl;
        }

        return matches;
    };

    virtual void compute(Mutation::Mixture& mix, Eigen::VectorXd& v) = 0;

private:
    double m_tolerance;
};

// Macro for adding a new comparison function
#define ADD_FUNCTION(__NAME__,__CODE__)\
class __NAME__ : public CompareFunc\
{\
public:\
    __NAME__(double tol) : CompareFunc(tol) {}\
protected:\
    void compute(Mutation::Mixture& mix, Eigen::VectorXd& v) {\
        __CODE__;\
    };\
};\
Mutation::Utilities::Config::ObjectProvider<\
    __NAME__, CompareFunc> op_##__NAME__(#__NAME__);

// Available comparison functions
ADD_FUNCTION(electron_thermal_conductivity_3,  v(0) = mix.electronThermalConductivity(3))
ADD_FUNCTION(electron_thermal_conductivity_2,  v(0) = mix.electronThermalConductivity(2))
ADD_FUNCTION(heavy_thermal_conductivity,     v(0) = mix.heavyThermalConductivity())
ADD_FUNCTION(internal_thermal_conductivity,  v(0) = mix.internalThermalConductivity(mix.T()))
ADD_FUNCTION(sigma_1st_order,                v(0) = mix.electricConductivity(1))
ADD_FUNCTION(sigma_2nd_order,                v(0) = mix.electricConductivity(2))
ADD_FUNCTION(heavy_thermal_diffusion_ratios, mix.heavyThermalDiffusionRatios(v.data()))
ADD_FUNCTION(viscosity,                      v(0) = mix.viscosity())

#endif // TESTS_COMPARE_FUNC_H
