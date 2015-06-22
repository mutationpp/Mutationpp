




#ifndef NUMERICS_LU_H
#define NUMERICS_LU_H

#include <cassert>
#include <vector>

#include "NumConst.h"
#include "Matrix.h"

namespace Mutation {
    namespace Numerics {

template <typename E>
LU(const MatExprt<Real, E>& A, bool do_nothing = false);

const Matrix<Real> &L();

const Matrix<Real> &U();



    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_LU_H
