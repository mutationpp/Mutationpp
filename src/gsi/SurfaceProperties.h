#ifndef SURFACEPROPERTIES_H
#define SURFACEPROPERTIES_H

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProperties{

public:
    typedef const Mutation::Utilities::IO::XmlElement& ARGS;

    SurfaceProperties( ARGS l_node_surf_props ){ }
    virtual ~SurfaceProperties(){ }

};

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // SURFACEPROPERTIES_H
