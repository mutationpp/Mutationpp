/**
 * @file Numerics.h
 *
 * @brief Convenience header to include all Numerics headers.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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

/**
 * @defgroup linsol Linear Solvers
 * Lists the available linear solvers in Mutation++ which solve the system
 * \f$Ax = b\f$ given the Numerics::Matrix \f$A\f$, and the vectors \f$x\f$ and
 * \f$b\f$.
 * @{
 * @defgroup itsol Iterative Solvers
 * Linear solvers based on iterative methods.
 *
 * @defgroup dirsol Direct Solvers
 * Linear solvers based on direct methods.
 * @}
 */

#ifndef NUMERICS_NUMERICS_H
#define NUMERICS_NUMERICS_H

#include "NumConst.h"
#include "Matrix.h"
#include "LeastSquares.h"
#include "lp.h"
#include "cg.h"
#include "gmres.h"
#include "LDLT.h"
#include "NewtonSolver.h"

namespace Mutation {
    namespace Numerics {


typedef Vector<double>           RealVector;
typedef VectorWrapper<double>    RealVecWrapper;
typedef Matrix<double>           RealMatrix;
typedef SymmetricMatrix<double>  RealSymMat;

typedef NumConst<double>         RealConsts;

    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_NUMERICS_H
