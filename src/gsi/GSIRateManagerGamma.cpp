#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerGamma : public GSIRateManager
{
public:

    GSIRateManagerGamma(DataGSIRateManager args)
        : GSIRateManager(args),
		  m_ns(args.s_thermo.nSpecies()),
		  m_nr(args.s_reactions.size()),
          v_wall_reaction_rate_constant(m_nr),
		  v_work(m_ns)
    {
        for (int i_reac = 0; i_reac < m_nr; ++i_reac){
            m_reactants.addReaction(i_reac, args.s_reactions[i_reac]->getReactants());
            m_irr_products.addReaction(i_reac, args.s_reactions[i_reac]->getProducts());
        }
    }

//=============================================================================

    ~GSIRateManagerGamma(){}

//=============================================================================

    /**
     * @toadddoc
     */
    Eigen::VectorXd computeRate()
    {
        // Get reaction rate constant
        for (int i_reac = 0; i_reac < m_nr; ++i_reac){
            v_wall_reaction_rate_constant(i_reac) =
            		v_reactions[i_reac]->getRateLaw()->forwardReactionRateCoefficient(
            				m_wall_state.getWallRhoi(), m_wall_state.getWallT());
        }

        // Constant rate times densities of species
        v_work.setZero();
        m_reactants.incrSpecies(v_wall_reaction_rate_constant, v_work);
        m_irr_products.decrSpecies(v_wall_reaction_rate_constant, v_work);

        // Multiply by molar mass
        return v_work.cwiseProduct(m_thermo.speciesMw().matrix());

    }

//=============================================================================

private:

    const int m_ns;
    const int m_nr;

    Eigen::VectorXd v_wall_reaction_rate_constant;
    Eigen::VectorXd v_work;

    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;

//=============================================================================

}; //class GSIRateManagerGamma

Mutation::Utilities::Config::ObjectProvider<GSIRateManagerGamma, GSIRateManager> gsi_rate_manager_gamma("gamma");
Mutation::Utilities::Config::ObjectProvider<GSIRateManagerGamma, GSIRateManager> gsi_rate_manager_gamma_energy("gamma_energy");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
