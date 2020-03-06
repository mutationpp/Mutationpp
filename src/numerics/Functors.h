/**
 * @file Functors.h
 *
 * @brief Defines several functors which are used throughout the project.
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

#ifndef NUMERICS_FUNCTORS
#define NUMERICS_FUNCTORS

#include <cmath>

namespace Mutation {
    namespace Numerics {


#define UNARY_FUNCTOR(__name__,__op__) \
template <typename T> \
class __name__ \
{ \
public: \
    inline T operator()(const T& x) const { return __op__; } \
};

/// Unary functor class f(x) = log(x).
UNARY_FUNCTOR(Log, std::log(x))

/// Unary functor class f(x) = exp(x).
UNARY_FUNCTOR(Exp, std::exp(x))

/// Unary functor class f(x) = sqrt(x).
UNARY_FUNCTOR(Sqrt, std::sqrt(x))

#undef UNARY_FUNCTOR

#define BINARY_FUNCTOR(__name__,__op__) \
template <typename T> \
class __name__ \
{ \
public: \
    inline T operator()(const T& x, const T& y) const { return __op__; } \
};

/// Binary functor class f(x,y) = x + y.
BINARY_FUNCTOR(Plus, x + y);

/// Binary functor class f(x,y) = x + abs(y).
BINARY_FUNCTOR(PlusAbs, x + std::abs(y))

/// Binary functor class f(x,y) = max(x, abs(y)).
BINARY_FUNCTOR(MaxAbs, std::max(x, std::abs(y)))

/// Binary functor class f(x,y) = max(x, y).
BINARY_FUNCTOR(Max, std::max(x, y))

/// Binary functor class f(x,y) = min(x, y).
BINARY_FUNCTOR(Min, std::min(x, y))

/// Binary functor class f(x,y) = x + y*y.
BINARY_FUNCTOR(PlusSquare, x + y*y);

#undef BINARY_FUNCTOR

#define EQUALITY_FUNCTOR(__name__,__op__) \
template <typename T> \
class __name__ \
{\
public:\
    inline void operator()(T& x, const T& y) const { x __op__ y; }\
};

/// Equality functor x = y
EQUALITY_FUNCTOR(Equals, =)

/// Equality functor x = x + y
EQUALITY_FUNCTOR(PlusEquals, +=)

/// Equality functor x = x - y
EQUALITY_FUNCTOR(MinusEquals, -=)

/// Equality functor x = x + a*y
template <typename T>
class PlusEqualsTimes
{
public:
    PlusEqualsTimes(const T& a)
        : m_a(a)
    { }
    inline void operator () (T& x, const T& y) const { x += m_a * y; }
private:
    const T& m_a;
};

/// Equality functor x = x - a*y
template <typename T>
class MinusEqualsTimes
{
public:
    MinusEqualsTimes(const T& a)
        : m_a(a)
    { }
    inline void operator () (T& x, const T& y) const { x -= m_a * y; }
private:
    const T& m_a;
};

/// Equality functor x = y / a
template <typename T>
class EqualsYDivAlpha
{
public:
    EqualsYDivAlpha(const T& a)
        : m_a(a)
    {}
    inline void operator()(T& x, const T& y) const { x = y / m_a; }
private:
    const T m_a;
};

/// Equality functor x += y / a
template <typename T>
class PlusEqualsYDivAlpha
{
public:
    PlusEqualsYDivAlpha(const T& a)
        : m_a(a)
    {}
    inline void operator()(T& x, const T& y) const { x += y / m_a; }
private:
    const T m_a;
};


    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_FUNCTORS
