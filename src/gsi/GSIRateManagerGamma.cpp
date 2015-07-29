#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerGamma : public GSIRateManager{

public:
    GSIRateManagerGamma( DataGSIRateManager l_data_rate_manager ) 
                    : GSIRateManager( l_data_rate_manager ),
                      v_wall_reaction_rate_constant( l_data_rate_manager.s_reactions.size(), 0.0 ) {

    for ( int i_reac = 0 ; i_reac < l_data_rate_manager.s_reactions.size() ; ++i_reac ){

        m_reactants.addReaction( i_reac, l_data_rate_manager.s_reactions[i_reac]->getReactants() );
        m_irr_products.addReaction( i_reac, l_data_rate_manager.s_reactions[i_reac]->getProducts() );
    
    }

}

//=============================================================================

    ~GSIRateManagerGamma(){ }

//=============================================================================

    void computeRate( Mutation::Numerics::RealVector& lv_mass_prod_rate ){ 

        // Get reaction rate constant
        for ( int i_reac = 0; i_reac < v_wall_reaction_rate_constant.size() ; ++i_reac ) {
            v_wall_reaction_rate_constant(i_reac) =  v_reactions[i_reac]->getRateLaw()->forwardReactionRate( m_surf_descr.getWallRhoi(), m_surf_descr.getWallT()); 
        }

        // Constant rate times densities of species
//        std::fill( lv_mass_prod_rate.begin(), lv_mass_prod_rate.end() , 0.0);
        for ( int i_sp = 0 ; i_sp < lv_mass_prod_rate.size() ; ++i_sp ){ lv_mass_prod_rate( i_sp ) = 0.E0; }
        m_reactants.incrSpecies( v_wall_reaction_rate_constant, lv_mass_prod_rate );
        m_irr_products.decrSpecies( v_wall_reaction_rate_constant, lv_mass_prod_rate );

        // Multiply by molar mass
        for ( int i_sp = 0 ; i_sp < lv_mass_prod_rate.size() ; ++i_sp ){
            lv_mass_prod_rate(i_sp) *= m_thermo.speciesMw(i_sp);
        }

    }

//=============================================================================
private:
    Mutation::Numerics::RealVector v_wall_reaction_rate_constant;

    GSIStoichiometryManagerGamma m_reactants;
    GSIStoichiometryManagerGamma m_irr_products;


//=============================================================================

};

Mutation::Utilities::Config::ObjectProvider<GSIRateManagerGamma, GSIRateManager> gsi_rate_manager_gamma("gamma");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
