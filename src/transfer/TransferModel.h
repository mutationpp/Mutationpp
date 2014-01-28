#ifndef TRANSFER_TRANSFER_MODEL_H
#define TRANSFER_TRANSFER_MODEL_H

#include "Thermodynamics.h"
#include "MillikanWhite.h"

namespace Mutation {
    namespace Transfer {

class TransferModel
{
public:
    virtual ~TransferModel() { }
    virtual void source(double* const p_source) = 0;
};

/**
 * Represents a coupling between vibrational and translational energy modes.
 */
class TransferModelVT : public TransferModel
{
public:

    TransferModelVT(const Thermodynamics::Thermodynamics& thermo)
        : m_mw(thermo)
    { }
    
    virtual void source(double* const p_source)
    {
    
    }

private:

    MillikanWhite m_mw;
};

    } // namespace Transfer
} // namespace Mutation

#endif // TRANSFER_TRANSFER_MODEL_H
