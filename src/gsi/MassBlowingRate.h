#ifndef MASSBLOWINGRATE_H
#define MASSBLOWINGRATE_H

#include "Utilities.h"

#include "DataMassBlowingRate.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class MassBlowingRate{
public:
    typedef const DataMassBlowingRate& ARGS;

    virtual ~MassBlowingRate(){ }

    virtual double computeBlowingFlux() = 0;
};


    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // MASSBLOWINGRATE_H
