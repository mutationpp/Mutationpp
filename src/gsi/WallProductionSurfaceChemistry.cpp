#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "WallProductionTerms.h"

#include "GSIRateManager.h"
#include "GSIReaction.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;

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
        DataGSIReaction data_gsi_reaction = { args.s_thermo,
                                              args.s_transport,
        		                              args.s_surf_props } ;
        for (data_gsi_reaction.s_iter_reaction = args.s_node_prod_terms.begin();
             data_gsi_reaction.s_iter_reaction != args.s_node_prod_terms.end();
             ++data_gsi_reaction.s_iter_reaction)
        {
            data_gsi_reaction.s_iter_reaction->
                getAttribute("type", m_reaction_type);

            addReaction(Factory<GSIReaction>::create(
                m_reaction_type, data_gsi_reaction));
        }

        // Assigning the mp_rate_manager pointer
        DataGSIRateManager data_rate_manager = { args.s_thermo,
                                                 args.s_surf_props,
                                                 args.s_wall_state,
                                                 mv_reaction };
        mp_rate_manager = Factory<GSIRateManager>::create(
            args.s_gsi_mechanism, data_rate_manager);
    }

//==============================================================================

    ~WallProductionSurfaceChemistry()
    {
        for (std::vector<GSIReaction*>::iterator iter = mv_reaction.begin();
        	 iter != mv_reaction.end(); ++iter){
        	delete (*iter);
        }
        mv_reaction.clear();

        if (mp_rate_manager != NULL){ delete mp_rate_manager; }
    }

//==============================================================================

    void productionRate(VectorXd& v_mass_prod_rate)
    {
        VectorXd v_wrk = mp_rate_manager->computeRate();

        v_mass_prod_rate.setZero();
        v_mass_prod_rate.head(m_ns) = v_wrk;
    }

//==============================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

private:
    const size_t m_ns;
    std::string m_tag;

    std::string m_reaction_type;

    std::vector<GSIReaction*> mv_reaction;
    GSIRateManager* mp_rate_manager;

    void addReaction(GSIReaction* l_gsi_reaction){
        mv_reaction.push_back(l_gsi_reaction);
    }

};

ObjectProvider<
    WallProductionSurfaceChemistry, WallProductionTerms>
    wall_production_surface_chemistry("surface_chemistry");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
