
#include "Mixture.h"
#include "StateModel.h"

namespace Mutation {

//==============================================================================

Mixture::Mixture(const MixtureOptions& options)
    : Thermodynamics::Thermodynamics(
        options.getSpeciesDescriptor(),
        options.getThermodynamicDatabase(),
        options.getStateModel()),
      Transport(
        *this,
        options.getViscosityAlgorithm(),
        options.getThermalConductivityAlgorithm(),
        options.loadTransport()),
      Kinetics(
        static_cast<const Thermodynamics&>(*this),
        options.getMechanism())
{
    // Set default composition if available
    if (options.hasDefaultComposition())
        setDefaultComposition(options.getDefaultComposition());
    
    // Instatiate a new energy transfer model
    mp_transfer = state()->createTransferModel(*this, *this, *this);
}

//==============================================================================

} // namespace Mutation
