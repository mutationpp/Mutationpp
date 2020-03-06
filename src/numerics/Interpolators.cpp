/**
 * @file Interpolators.cpp
 *
 * @brief Registers the concrete Interpolator types to support self
 * registration.
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

#include "Interpolators.h"
#include "Utilities.h"

namespace Mutation {
    namespace Numerics {

// ChebyshevInterpolator
Utilities::Config::ObjectProvider<
    ChebyshevInterpolator<double>, Interpolator<double> > chebyshev_d("Chebyshev");
Utilities::Config::ObjectProvider<
    ChebyshevInterpolator<float>,  Interpolator<float> >  chebyshev_f("Chebyshev");

// LinearInterpolator
Utilities::Config::ObjectProvider<
    LinearInterpolator<double>, Interpolator<double> > linear_d("Linear");
Utilities::Config::ObjectProvider<
    LinearInterpolator<float>,  Interpolator<float> >  linear_f("Linear");

// MCHInterpolator
Utilities::Config::ObjectProvider<
    MCHInterpolator<double>, Interpolator<double> > mch_d("MonotoneCubic");
Utilities::Config::ObjectProvider<
    MCHInterpolator<float>, Interpolator<float> >   mch_f("MonotoneCubic");

    } // Numerics
} // Mutation
