#ifndef NUMERICS_NUMCONST_H
#define NUMERICS_NUMCONST_H

#include <limits>

namespace Mutation {
    namespace Numerics {


/**
 * Stores numeric constants for a template type which are used throughout the
 * algorithms in the numerics directory.
 */
template <typename T>
class NumConst
{
public:
    static const T zero;
    static const T one;
    static const T two;
    static const T eps; ///< Machine round-off
};

template <typename T>
const T NumConst<T>::zero = static_cast<T>(0);

template <typename T>
const T NumConst<T>::one  = static_cast<T>(1);

template <typename T>
const T NumConst<T>::two  = static_cast<T>(2);

template <typename T>
const T NumConst<T>::eps  = std::numeric_limits<T>::epsilon();


    } // namespace Numerics
} // namespace Mutation

#endif // NUMERICS_NUMCONST_H
