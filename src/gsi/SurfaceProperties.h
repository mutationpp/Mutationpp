#ifndef SURFACE_PROPERTIES_H
#define SURFACE_PROPERTIES_H

#include "DataSurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProperties{

public:
	/**
	 *
	 */
    typedef const DataSurfaceProperties& ARGS;

    /**
     *
     */
    SurfaceProperties(ARGS args){ }

    /**
     *
     */
    virtual ~SurfaceProperties(){ }

public:
    // FRC DATA
    virtual int speciesIndexWall( const std::string& str_sp ) const = 0;
    virtual int nSpeciesWall() const = 0;

    virtual int nSites() const = 0;
    virtual double nTotalSites() const = 0;

    virtual double fracSite( const int& i_site ) const = 0;
    virtual int nSpeciesSite( const int& i_site ) const = 0;

    // EQUILIBRIUM DATA @BD Improve names?
    virtual bool surfaceConstraint() const {return false;}
    virtual std::string surfaceSpecies() const {return "";}
    virtual double BprimePyro() const {return 0;}
    virtual void getSurfaceCompositions(
                                     Eigen::VectorXd & lv_pyrolisysComposition, 
                                     Eigen::VectorXd & lv_charComposition){}

}; // class SurfaceProperties

    } // namespace GasSurfaceInteraction 
} // namespace Mutation

#endif // SURFACE_PROPERTIES_H
