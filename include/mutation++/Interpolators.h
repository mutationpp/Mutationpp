/**
 * @file Interpolators.h
 *
 * @brief Provides Interpolator base class and derived interpolator types.
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

#ifndef INTERPOLATORS_H
#define INTERPOLATORS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace Mutation {
    namespace Numerics {

/// Self registering base class for all Interpolators.
template <typename T>
class Interpolator
{
public:

    /// Returns name of this type.
    static std::string typeName() { return "Interpolator"; }

    // Arguments for self registering constructor
    struct ARGS {
        ARGS(
            const T* const arg1, const T* const arg2, int arg3) :
                x(arg1), y(arg2), n(arg3) {}
        const T* const x;
        const T* const y;
        int n;
    };

    /// Constructor used in self registration.
    Interpolator(ARGS args) { }

    /// Destructor
    virtual ~Interpolator() { }

    /// Returns interpolated function.
    virtual T operator() (const T& x) = 0;
}; // Interpolator

//==============================================================================

/**
 * Chebyshev polynomial interpolator.
 *
 * Authors: J.B. Scoggins, translated from Matlab code of Fernando Miro Miro
 */
template <typename T>
class ChebyshevInterpolator : public Interpolator<T>
{
public:

    typedef Eigen::Array<T,-1,1> ArrayType;

    /// Takes a list of points and generates the interpolation function.
    ChebyshevInterpolator(typename Interpolator<T>::ARGS args) :
        Interpolator<T>(args),
        m_points(args.n-1), m_eta_map(args.n-1), m_y_map(args.n-1),
        m_phi(args.n-1)
    {
        const T* const x = args.x;
        const T* const y = args.y;
        const int n = args.n;
        const double PI = std::acos(-1.0);

        // Minimum, median, and maximum
        m_min = x[0];
        m_mid = x[(n-1)/2];
        m_max = x[n-1];

        // Compute the Chebyshev-Gauss-Lobatto grid
        ArrayType cgl(m_points);
        for (int i = 0; i < m_points; ++i)
            cgl(i) = std::cos(PI*i/(m_points-1.0));

        // Compute the new
        T mid = m_mid - m_min;
        T max = m_max - m_min;
        T a   = mid*max/(max - 2.0*mid);
        T b   = 1.0 + 2.0*a/max;
        ArrayType x_map = a*(1.0 + cgl)/(b - cgl) + m_min;

        m_eta_map = (m_mid - x_map)/(2.0*m_mid*x_map/m_max - x_map - m_mid);

        // Interpolating points onto new map
        int k;
        for (int i = 0; i < m_points; ++i) {
            k = 1;
            while (x[k] < x_map(i) && k < m_points) ++k;
            m_y_map(i) = (x_map(i) - x[k-1])*(y[k] - y[k-1])/(x[k] - x[k-1]) + y[k-1];
        }
    }

    int nPoints() const { return m_points; }

    /// Interpolates the function at x.
    T operator() (const T& x)
    {
        T eta = (m_mid - x)/(2.0*m_mid*x/m_max - x - m_mid);

        for (int j = 0; j < m_points; ++j) {
            m_phi(j) = 1.0;
            for (int k = 0; k < j; ++k)
                m_phi(j) *= (eta - m_eta_map(k))/(m_eta_map(j) - m_eta_map(k));
            for (int k = j+1; k < m_points; ++k)
                m_phi(j) *= (eta - m_eta_map(k))/(m_eta_map(j) - m_eta_map(k));
        }

        return (m_phi*m_y_map).sum();
    }

private:

    int m_points;

    ArrayType m_eta_map;
    ArrayType m_y_map;
    ArrayType m_phi;

    T m_min;
    T m_mid;
    T m_max;
}; // ChebyshevInterpolator

//==============================================================================

/// Linear interpolator.
template <typename T>
class LinearInterpolator : public Interpolator<T>
{
public:

    typedef Eigen::Array<T,-1,1> ArrayType;

    LinearInterpolator(typename Interpolator<T>::ARGS args) :
        Interpolator<T>(args),
        m_x(Eigen::Map<const ArrayType>(args.x, args.n)),
        m_y(Eigen::Map<const ArrayType>(args.y, args.n))
    { }

    int nPoints() const { return m_x.size(); }

    /// Interpolates the function at x.
    T operator() (const T& x)
    {
        int i = 1;
        while (m_x[i] < x && i < m_x.size()-1) i++;

        // Interpolate
        return (m_y[i]-m_y[i-1])*(x-m_x[i])/(m_x[i]-m_x[i-1])+m_y[i];
    }

private:

    const ArrayType m_x;
    const ArrayType m_y;

}; // LinearInterpolator

//==============================================================================

/// Monotone Cubic Hermite interpolator.
template <typename T>
class MCHInterpolator : public Interpolator<T>
{
public:
    /// Takes a list of points and generates the interpolation function.
    MCHInterpolator(typename Interpolator<T>::ARGS args) :
        Interpolator<T>(args),
        m_points(args.n), m_x(args.n), m_y(args.n), m_c1(args.n), m_c2(args.n),
        m_c3(args.n)
    {
        const T* const x = args.x;
        const T* const y = args.y;
        const int n = args.n;

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
    T operator() (const T& x)
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
}; // MCHInterpolator

//==============================================================================

    } // Numerics
} // Mutation

#endif // INTERPOLATORS_H
