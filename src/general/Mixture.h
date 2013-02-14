#ifndef MUTATION_MIXTURE_H
#define MUTATION_MIXTURE_H

#include "Thermodynamics.h"
#include "Kinetics.h"
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
    : public Thermodynamics::Thermodynamics, public Transport, public Kinetics
{
public:

    /**
     * Constructs a Mixture using the options given.  Note that a mixture file
     * name can also be given instead and the options will be loaded 
     * accordingly.
     *
     * @see MixtureOptions
     */
    Mixture(const MixtureOptions& options)
        : Thermodynamics::Thermodynamics(
             options.getSpeciesNames(), 
             options.getThermodynamicDatabase(),
             options.getStateModel()),
          Transport(
             *this, 
             options.getViscosityAlgorithm(),
             options.getThermalConductivityAlgorithm()),
          Kinetics(
             static_cast<const Thermodynamics&>(*this),
             options.getMechanism())
    { 
        if (options.hasDefaultComposition())
            setDefaultComposition(options.getDefaultComposition());
    }
    
    /** 
     * Destructor.
     */
    ~Mixture() {}
    
}; // class Mixture

} // namespace Mutation

#endif // MUTATION_MIXTURE_H

