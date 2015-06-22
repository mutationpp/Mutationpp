#ifndef BUILDER_H
#define BUILDER_H

#include "Thermodynamics.h"
#include "CatalysisRateManager.h"

namespace Mutation {
    namespace gsi {

class BuilderGasSurfaceInteraction {
public:

    BuilderGasSurfaceInteraction( const std::string& l_catalytic_model, const Mutation::Thermodynamics::Thermodynamics& thermo );
    ~BuilderGasSurfaceInteraction( );

//        CatalysisRateManager* getCatalysisRateManagerModel( CatalysisSurfaceProperties* const lp_surface_properties, const std::vector<CatalysisReaction>& l_catalytic_reactions );
//        CatalysisSurfaceProperties* getCatalysisSurfacePropertiesModel( const std::string& l_surface_file );
//        WallState* getWallStateModel( CatalysisSurfaceProperties* const lp_surface_properties );
//    //    FactoryRateLaws* getFactoryRateLaws
//    
//    private:
//        void giveErrorMessageandExit();
//    
//        enum string_catalytic_model {
//            gamma,
//            finite_rate_chemistry,
//            error
//        };
//    
//        string_catalytic_model convertStringtoEnumCatalyticModel (const std::string& l_catalytic_model);
//    
//    
//    private:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    std::string m_catalytic_model;

};

    } // namespace gsi
} // namespace Mutation

#endif // BUILDER_H
