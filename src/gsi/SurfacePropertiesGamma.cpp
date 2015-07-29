#include "AutoRegistration.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfacePropertiesGamma : public SurfaceProperties {

public:
    SurfacePropertiesGamma( const Mutation::Utilities::IO::XmlElement& l_node_surf_props ) : SurfaceProperties( l_node_surf_props ) { }
    ~SurfacePropertiesGamma(){ }

};

Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesGamma, SurfaceProperties> surface_properties_gamma("gamma");


    } // namespace GasSurfaceInteraction
} // namespace Mutation 
