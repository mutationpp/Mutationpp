
#include "ThermoDB.h"
#include "Species.h"

#include <algorithm>

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ThermoDB::ThermoDB(ThermoDB::ARGS species)
    : m_ns(species.size())
{ }

//==============================================================================

// Used in ThermoDB::cv
struct MinusOne {
    void operator () (double& x) const { x -= 1.0; }
} MinusOne;

void ThermoDB::cv(
    double Th, double Te, double Tr, double Tv, double Tel, 
    double* const cv = NULL, double* const cvt = NULL, double* const cvr = NULL,
    double* const cvv = NULL, double* const cvel = NULL)
{
    // Compute Cp/Ru
    cp(Th, Te, Tr, Tv, Tel, cv, cvt, cvr, cvv, cvel);

    // Cv/Ru = Cp/Ru - 1
    if (cv   != NULL) std::for_each(  cv+0,   cv+m_ns, MinusOne);        
    if (cvt  != NULL) std::for_each( cvt+0,  cvt+m_ns, MinusOne);
    if (cvr  != NULL) std::for_each( cvr+0,  cvr+m_ns, MinusOne);
    if (cvv  != NULL) std::for_each( cvv+0,  cvv+m_ns, MinusOne);
    if (cvel != NULL) std::for_each(cvel+0, cvel+m_ns, MinusOne);
}

//==============================================================================

} // namespace Thermodynamics
} // namespace Mutation

