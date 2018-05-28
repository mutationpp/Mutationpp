#include "AutoRegistration.h"
#include "Utilities.h"

#include "DataGSIRateManager.h"
#include "DataGSIReaction.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionSurfaceChemistry : public WallProductionTerms
{
public:
    WallProductionSurfaceChemistry(ARGS args)
        : WallProductionTerms(args),
          m_ns(args.s_thermo.nSpecies()),
          mp_rate_manager(NULL),
          m_tag("surface_chemistry")
    {

        // Filling in the reaction vector
        DataGSIReaction l_data_gsi_reaction = { args.s_thermo,
                                                args.s_transport,
        		                                args.s_surf_props } ;
        for (l_data_gsi_reaction.s_iter_reaction = args.s_node_prod_terms.begin();
             l_data_gsi_reaction.s_iter_reaction != args.s_node_prod_terms.end();
             ++l_data_gsi_reaction.s_iter_reaction)
        {
            l_data_gsi_reaction.s_iter_reaction->getAttribute("type", m_reaction_type); /** @todo Give an error! */
            addReaction( Mutation::Utilities::Config::Factory<GSIReaction>::create(
                m_reaction_type, l_data_gsi_reaction));
        }

        // Assigning the mp_rate_manager pointer
        DataGSIRateManager l_data_rate_manager = { args.s_thermo,
                                                   args.s_surf_props,
												   args.s_wall_state,
                                                   v_reaction };
        mp_rate_manager = Mutation::Utilities::Config::Factory<GSIRateManager>::create(
            args.s_gsi_mechanism, l_data_rate_manager);

    }

//======================================================================================

    ~WallProductionSurfaceChemistry()
    {
        for (std::vector<GSIReaction*>::iterator iter = v_reaction.begin(); iter != v_reaction.end(); ++iter){
        	delete (*iter);
        }
        v_reaction.clear();

        if (mp_rate_manager != NULL){ delete mp_rate_manager; }
    }

//======================================================================================

    void productionRate(Eigen::VectorXd& lv_mass_prod_rate)
    {
        Eigen::VectorXd lv_wrk = mp_rate_manager->computeRate();

        lv_mass_prod_rate.setZero();
        lv_mass_prod_rate.head(m_ns) = lv_wrk;
    }

//======================================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

//======================================================================================

private:
    const size_t m_ns;
    std::string m_tag;

    std::string m_reaction_type;

    std::vector<GSIReaction*> v_reaction;
    GSIRateManager* mp_rate_manager;

    void addReaction(GSIReaction* l_gsi_reaction){
        v_reaction.push_back(l_gsi_reaction);
    }

}; // class WallProductionSurfaceChemistry

//======================================================================================

Mutation::Utilities::Config::ObjectProvider<WallProductionSurfaceChemistry, WallProductionTerms> wall_production_surface_chemistry("surface_chemistry");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
