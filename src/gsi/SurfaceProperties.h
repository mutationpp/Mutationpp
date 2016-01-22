#ifndef SURFACEPROPERTIES_H
#define SURFACEPROPERTIES_H

#include "DataSurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProperties{

public:
    typedef const DataSurfaceProperties& ARGS;

    SurfaceProperties( ARGS l_data_surf_props ){ }
    virtual ~SurfaceProperties(){ }

public:
    virtual int speciesIndexWall( const std::string& str_sp ) const = 0;
    virtual int nSpeciesWall() const = 0;

    virtual int nSites() const = 0;
    virtual double nTotalSites() const = 0;

    virtual double fracSite( const int& i_site ) const = 0;
    virtual int nSpeciesSite( const int& i_site ) const = 0;

};

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // SURFACEPROPERTIES_H
