#ifndef SURFACEPRODUCTIONTERMS_H
#define SURFACEPRODUCTIONTERMS_H

#include "GSIReaction.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProductionTerms {

public:
    SurfaceProductionTerms(){ }
    virtual ~SurfaceProductionTerms(){ }

    virtual void productionRate() = 0;

};

//=========================================================

class SurfaceProductionRate : public SurfaceProductionTerms {

public:
    SurfaceProductionRate(){ }
    virtual ~SurfaceProductionRate(){ }

    void productionRate(){ };

private:
    std::vector<GSIReaction> v_reaction;

};

//=========================================================

class BulkProductionRate : public SurfaceProductionTerms {

public:
    BulkProductionRate(){ }
    virtual ~BulkProductionRate(){ }

    void productionRate(){ };

};

//=========================================================

class FailureProductionRate : public SurfaceProductionTerms {

public:
    FailureProductionRate(){ }
    virtual ~FailureProductionRate(){ }

    void productionRate(){ };

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation



#endif // SURFACEPRODUCTIONTERMS_H
