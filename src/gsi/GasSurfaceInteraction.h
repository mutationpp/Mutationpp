#ifndef GSI_H
#define GSI_H

#include "Utilities.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "GSIReaction.h"
#include "CatalysisRateManager.h"
#include "SurfaceProperties.h"
#include "WallSolver.h"
//#include "BuilderGasSurfaceInteraction.h"

namespace Mutation {
    namespace gsi {
      
class GasSurfaceInteraction {
public:

    /**
     * Constructor which takes a constant reference to a Thermodynamics object,
     * a reference to a Transport object for the computation of the mass balance
     * at the wall and a gsi mechanism file.
     */
    GasSurfaceInteraction(Mutation::Thermodynamics::Thermodynamics& thermo,
                          Mutation::Transport::Transport& transport,
                          const std::string gsi_catalysis_mechanism_file);

    /**
     * Destructor.
     */
    ~GasSurfaceInteraction();
    
    /**
     * Returns the number of reactions in the mechanism
     */
    size_t nCatalyticReactions() const {return v_catalytic_reactions.size();}

    /**
     * Returns a pointer containing the gamma coefficients for reaction i in the mechanism
     */
//    void getGammaCoefficientforReactioni(const int& i_reaction, double* p_gamma_coefficients);

    /**
     * Computes the chemical source terms for the gas phase species that are produced according to the
     * chosen catalytic model.
     * In order to use this model the wall state should have been set by calling the setWallState function.
     * The signs for the production rates are from the wall frame of reference, meaning positive for the
     * species impinging to the wall and negative for the species leaving the wall.
     *
     * @param p_wdot - on return, the species production rates in kg/m^3-s
     *
     */
    void netGSIProductionRates(double * const p_wdot);

    /** @todo
     * Solves for composition of the gas species on the surface by solving the mass balance equation:
     * \f$ J_i(\rho_i, \mathbf{n} \cdot \nabla) \f$
     *
     * Overloaded with the order of computing the gradient at the wall.
     */
    void solveWallMassBalance(double* const p_rhoi, const double* const p_rhoi1, 
                              const double* const p_rhoie, const double& dx_length_mole_fraction_grad);

    /**
     * Call in order to set the state at the wall. Necessary to be called before asking for the netGSIProductionRates;
     * @param p_rhoi partial densities at the in kg/m^3
     * @param p_rhoie temperatures of the wall. Only the first one, translational, is taken currently into account in K
     *
     * @todo Update setWallState function using state models.
     */

    void setWallState(const double* const p_rhoi, const double* const p_rhoie);

private:
    /**
     * Adds a new reaction to the GasSurfaceInteraction object.
     */
    void addCatalyticReaction(const CatalysisReaction& catalytic_reaction);

    /**
     * Closes reactions vector. The control for mass and elements conservation is performed here.
     * Trying to addCatalyticReaction() after the closeGSIReactions() is called results to an error.
     *
     * @todo To be implemented
     */
    void closeGSIReactions(const bool validate_gsi_mechanism);

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Transport::Transport& m_transport;

    std::vector<CatalysisReaction> v_catalytic_reactions;
    std::string m_category;
    std::string m_catalytic_model;
    std::string m_surface;

//    BuilderGasSurfaceInteraction* mp_builder_gassurfaceinteraction;
    CatalysisRateManager* mp_catalysis_rates;
    CatalysisSurfaceProperties* mp_surface_properties;
    WallState* mp_wall_state;
    WallSolver* mp_wall_solver;

    bool m_wall_state_set;
};

    } // namespace gsi
} // namespace Mutation

#endif // GSI_H
