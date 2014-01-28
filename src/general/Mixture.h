#ifndef MUTATION_MIXTURE_H
#define MUTATION_MIXTURE_H

#include "Thermodynamics.h"
#include "Kinetics.h"
#include "TransferModel.h"
#include "Transport.h"
#include "MixtureOptions.h"

namespace Mutation {

/**
 * Packages all of the Mutation++ functionality pertaining to a mixture into a
 * single class.  A Mixture object simply inherits all functionalities from the
 * Thermodynamics, Transport, and Kinetics classes.  A Mixture object can be
 * constructed using a mixture file name or a MixtureOptions object.
 *
 * @see Thermodynamics
 * @see Transport
 * @see Kinetics
 */
class Mixture
    : public Thermodynamics::Thermodynamics, 
      public Transport::Transport, 
      public Kinetics::Kinetics
{
public:

    /**
     * Constructs a Mixture using the options given.  Note that a mixture file
     * name can also be given instead and the options will be loaded 
     * accordingly.
     *
     * @see MixtureOptions
     */
    Mixture(const MixtureOptions& options);
    
    /** 
     * Destructor.
     */
    ~Mixture()
    {
        if (mp_transfer != NULL) delete mp_transfer;
    }
    
    /**
     * Provides energy transfer source terms based on the current state of the
     * mixture.
     */
    void energyTransferSource(double* const p_source) {
        if (mp_transfer != NULL)
            mp_transfer->source(p_source);
    }

private:

    Transfer::TransferModel* mp_transfer;
    
}; // class Mixture

} // namespace Mutation

#endif // MUTATION_MIXTURE_H

