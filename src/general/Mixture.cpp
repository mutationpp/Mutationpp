
#include "Mixture.h"
#include "StateModel.h"

/**
 * @todo Check if this is needed
 */

#include "GasSurfaceInteraction.h"

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
        options.getMechanism()),
      GasSurfaceInteraction(
        static_cast<const Thermodynamics&>(*this),
        options.getGSIMechanism() 
        /** @todo ablation , options.get*/)
{
    // Set default composition if available
    if (options.hasDefaultComposition())
        setDefaultComposition(options.getDefaultComposition());
    
    // Instatiate a new energy transfer model
    mp_transfer = state()->createTransferModel(*this, *this, *this);
    
    //Initializing gsi module
    /**
     * @todo remove GSIFLAG. Make it compile by default 
     */ 
//    #ifdef GSIFLAG
//        gsi::GasSurfaceInteraction* mp_gsi = new gsi::GasSurfaceInteraction(static_cast<const Thermodynamics&>(*this), options.getGSIMechanism() /** @todo ablation , options.get*/);
//    #endif
}

//==============================================================================

} // namespace Mutation
