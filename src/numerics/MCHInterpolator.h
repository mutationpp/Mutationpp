/**
 * @file MCHInterpolator.h
 *
 * @brief Provides a simple monotone cubic hermite interpolator object.
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

#ifndef MCH_INTERPOLATOR_H
#define MCH_INTERPOLATOR_H

#include <algorithm>
#include <vector>

namespace Mutation {
    namespace Numerics {

/// Monotone Cubic Hermite interpolator.
template <typename T>
class MCHInterpolator
{
public:
    /// Takes a list of points and generates the interpolation function.
    MCHInterpolator(const T* const x, const T* const y, int n) :
        m_points(n), m_x(n), m_y(n), m_c1(n), m_c2(n), m_c3(n)
    {
        // Store the x and y values
        for (int i = 0; i < n; ++i) {
            m_x[i] = x[i];
            m_y[i] = y[i];
        }

        // Get consecutive differences and slopes
        std::vector<T> dx(n);
        std::vector<T> dy(n);
        std::vector<T> m(n);
        for (int i = 0; i < n-1; ++i) {
            dx[i] = x[i+1] - x[i];
            dy[i] = y[i+1] - y[i];
            m[i]  = dy[i] / dx[i];
        }

        // Get degree-1 coefficients
        m_c1[0] = m[0];
        for (int i = 1; i < n-1; ++i) {
            if (m[i-1]*m[i] <= 0)
                m_c1[i] = T(0);
            else {
                T c = dx[i-1] + dx[i];
                m_c1[i] = 3*c/((c + dx[i])/m[i-1] + (c + dx[i-1])/m[i]);
            }
        }
        m_c1[n-1] = m[n-1];

        // Get degree-2 and degree-3 coefficients
        for (int i = 0; i < n-1; ++i) {
            T c = m_c1[i] + m_c1[i+1] - 2*m[i];
            m_c2[i] = (m[i] - m_c1[i] - c)/dx[i];
            m_c3[i] = c/(dx[i]*dx[i]);
        }
    }

    int nPoints() const { return m_points; }

    /// Interpolates the function at x.
    T operator() (const T& x) const
    {
        if (x >= m_x[m_points-1])
            return m_y[m_points-1];

        if (x <= m_x[0])
            return m_y[0];

        int i = static_cast<int>(
            std::lower_bound(&m_x[0], &m_x[0] + m_points - 1, x) - &m_x[0]) - 1;

        T d = x - m_x[i];
        return m_y[i] + ((m_c3[i]*d + m_c2[i])*d + m_c1[i])*d;
    }

private:

    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

private:

    int m_points;
    std::vector<T> m_x;
    std::vector<T> m_y;
    std::vector<T> m_c1;
    std::vector<T> m_c2;
    std::vector<T> m_c3;
};

    } // Numerics
} // Mutation

#endif // MCH_INTERPOLATOR_H
