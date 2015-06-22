#include <string>
#include <vector>

#include "BuilderGasSurfaceInteraction.h"

namespace Mutation {
    namespace gsi {

BuilderGasSurfaceInteraction::BuilderGasSurfaceInteraction( const std::string& l_catalytic_model, const Mutation::Thermodynamics::Thermodynamics& thermo )
                                                          : m_catalytic_model( l_catalytic_model ),
                                                            m_thermo( thermo )
{



}

BuilderGasSurfaceInteraction::~BuilderGasSurfaceInteraction( ){ }

//===============================================================================

//    CatalysisRateManager* BuilderGasSurfaceInteraction::getCatalysisRateManagerModel( CatalysisSurfaceProperties* const lp_surface_properties, const std::vector<CatalysisReaction>& l_catalytic_reactions ){
//    
//        switch ( convertStringtoEnumCatalyticModel(m_catalytic_model) )
//        {
//            case gamma:
//                return new CatalysisGammaRateManager( m_thermo, l_catalytic_reactions ) ;
//            case finite_rate_chemistry:
//                return new CatalysisFRCRateManager ( m_thermo, lp_surface_properties, l_catalytic_reactions );
//            default:
//                giveErrorMessageandExit();
//                return NULL;
//        }
//    
//    }
//    
//    //===============================================================================
//    
//    CatalysisSurfaceProperties* BuilderGasSurfaceInteraction::getCatalysisSurfacePropertiesModel( const std::string& l_surface_file ){
//    
//        switch ( convertStringtoEnumCatalyticModel( m_catalytic_model ) )
//        {
//            case gamma:
//                return NULL;
//            case finite_rate_chemistry:
//                return new CatalysisSurfaceProperties( l_surface_file, m_thermo );
//            default:
//                giveErrorMessageandExit();
//                return NULL;
//        }
//    
//    }
//    
//    //===============================================================================
//    
//    WallState* BuilderGasSurfaceInteraction::getWallStateModel( CatalysisSurfaceProperties* const lp_surface_properties ){
//    
//        switch ( convertStringtoEnumCatalyticModel(m_catalytic_model) )
//        {
//            case gamma:
//                return new WallState( m_thermo );
//            case finite_rate_chemistry:
//                return new WallStateFRC( *lp_surface_properties, m_thermo );
//            default:
//                giveErrorMessageandExit();
//                return NULL;
//        }
//    
//    }
//    
//    //===============================================================================
//    
//    void BuilderGasSurfaceInteraction::giveErrorMessageandExit(){
//    
//        std::cout << "The catalytic model " << m_catalytic_model << " has not been implemented yet!";
//        exit(1);
//    
//    }
//    
//    //===============================================================================
//    
//    BuilderGasSurfaceInteraction::string_catalytic_model BuilderGasSurfaceInteraction::convertStringtoEnumCatalyticModel ( const std::string& l_catalytic_model ){
//    
//        if ( l_catalytic_model == "gamma" ) return gamma;
//        if ( l_catalytic_model == "finite_rate_chemistry" ) return finite_rate_chemistry;
//    
//        giveErrorMessageandExit();
//        return error;
//    
//    }

    } // namespace gsi
} // namespace Mutation

