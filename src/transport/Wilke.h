/**
 * @file Wilke.h
 *
 * @brief Provides Wilke class.
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

#ifndef TRANSPORT_WILKE_H
#define TRANSPORT_WILKE_H

#include <Eigen/Dense>

namespace Mutation {
    namespace Transport {

/**
 * Base class for classes that compute transport properties using the Wilke
 * formulas.
 */
class Wilke
{
public:

    Wilke() { }

protected:

    template <typename E1, typename E2, typename E3>
    double wilke(
        const Eigen::ArrayBase<E1>& vals, const Eigen::ArrayBase<E2>& mass,
        const Eigen::ArrayBase<E3>& x)
    {
        const int ns = vals.size();
        double average = 0.0, sum, ratio, temp;

        for (int i = 0; i < ns; ++i) {
            sum = 0.0;

            for (int j = 0; j < ns; ++j) {
                if (i == j)
                    sum += x(j);
                else {
                    ratio = mass(i) / mass(j);
                    temp = 1.0 + std::sqrt(vals(i) / vals(j) / std::sqrt(ratio));
                    sum += x(j) * temp * temp / std::sqrt(8.0 * (1.0 + ratio));
                }
            }

            average += x(i) * vals(i) / sum;
        }

        return average;
    }

}; // class Wilke


    } // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_WILKE_H
