#include "AutoRegistration.h"
#include "Utilities.h"

#include "DataGSIRateManager.h"
#include "DataGSIReaction.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProductionRate : public WallProductionTerms {

//======================================================================================

public:
    SurfaceProductionRate( ARGS l_data_wall_prod_terms )
                         : WallProductionTerms( l_data_wall_prod_terms ),
                           mp_rate_manager( NULL ) {

        // Filling in the reaction vector
        DataGSIReaction l_data_gsi_reaction = { l_data_wall_prod_terms.s_thermo, l_data_wall_prod_terms.s_surf_props } ;
        for (l_data_gsi_reaction.s_iter_reaction = l_data_wall_prod_terms.s_node_prod_terms.begin() ; l_data_gsi_reaction.s_iter_reaction != l_data_wall_prod_terms.s_node_prod_terms.end() ; ++l_data_gsi_reaction.s_iter_reaction ) {
            l_data_gsi_reaction.s_iter_reaction->getAttribute("type", m_reaction_type ); /** @todo Give an error! */ 
            addReaction( Mutation::Utilities::Config::Factory<GSIReaction>::create( m_reaction_type, l_data_gsi_reaction ) );
        }

        // Assigning the mp_rate_manager pointer
        DataGSIRateManager l_data_rate_manager = { l_data_wall_prod_terms.s_thermo, l_data_wall_prod_terms.s_surf_props, l_data_wall_prod_terms.s_wall_state, v_reaction };
        mp_rate_manager = Mutation::Utilities::Config::Factory<GSIRateManager>::create( l_data_wall_prod_terms.s_gsi_mechanism, l_data_rate_manager );

    }

//======================================================================================

    ~SurfaceProductionRate(){

        for ( std::vector<GSIReaction*>::iterator iter = v_reaction.begin() ; iter != v_reaction.end() ; ++iter ) delete (*iter);
        
        v_reaction.clear();
        if ( mp_rate_manager != NULL ){ delete mp_rate_manager; }

    }

//======================================================================================

    void productionRate( Eigen::VectorXd& lv_mass_prod_rate ){ 

        mp_rate_manager->computeRate( lv_mass_prod_rate );

    }

//======================================================================================

private:
    std::string m_reaction_type;

    std::vector<GSIReaction*> v_reaction; /** @todo auto_ptr */
    GSIRateManager* mp_rate_manager;

    void addReaction( GSIReaction* l_gsi_reaction ){
        v_reaction.push_back( l_gsi_reaction );
    }

};

//======================================================================================

Mutation::Utilities::Config::ObjectProvider<SurfaceProductionRate, WallProductionTerms> surface_production_rate("surface_chemistry");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
