#ifndef CATRATEMANAGER_H
#define CATRATEMANAGER_H

#include "GSIReaction.h"
#include "GSIStoichiometryManager.h"

#include "AutoRegistration.h"

/**
 * This class organizes the GSIReactions and computes the necessary matrices and 
 * laws in order to compute in the end the rate law etc. 
 * @todo Fix this description.
 */

namespace Mutation{
    namespace gsi{
      
//=============================================================================
//=============================================================================
//=============================================================================

class CatalysisRateManager{

public:
    virtual ~CatalysisRateManager(){ }
    virtual CatalysisRateManager* clone() const = 0;

    virtual void computeRate(const WallState& r_wall_state, double* const p_prod) = 0;
};

//=============================================================================
//=============================================================================
//=============================================================================

class CatalysisRateManagerGamma : public CatalysisRateManager{

public:
    CatalysisRateManagerGamma( const Mutation::Thermodynamics::Thermodynamics& thermo, const std::vector<CatalysisReaction>& catalytic_reactions );
    ~CatalysisRateManagerGamma() { }

    CatalysisRateManagerGamma* clone() const {
        return new CatalysisRateManagerGamma(*this);
    }
    
    void computeRate(const WallState& r_wall_state, double* const p_prod);
    
private:
    const size_t m_ns;
    const size_t m_nr;

    Mutation::Numerics::RealVector v_wall_reaction_rate_constant;

    const std::vector<CatalysisReaction>& v_reactions;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    GSIStoichiometryManager_GammaModel v_reactants;
    GSIStoichiometryManager_GammaModel m_irr_products;
};

//=============================================================================
//=============================================================================
//=============================================================================

class CatalysisRateManagerFiniteRateChemistry : public CatalysisRateManager{

public:
    CatalysisRateManagerFiniteRateChemistry(const Mutation::Thermodynamics::Thermodynamics& thermo, const Mutation::gsi::CatalysisSurfaceProperties* mp_surf_props, const std::vector<Mutation::gsi::CatalysisReaction>& m_catalytic_reactions);
    ~CatalysisRateManagerFiniteRateChemistry();

    CatalysisRateManagerFiniteRateChemistry* clone() const{
        return new CatalysisRateManagerFiniteRateChemistry(*this);
    }

    void computeRate(const WallState& l_wall_state, double* const p_prod);
    
private:
    void solveSteadyStateCoverage(const WallState& r_wall_state);
    
    const size_t m_ns;
    const size_t m_nr;
    const size_t m_nsites;
    const size_t m_total_species_in_sites;

    Mutation::Numerics::RealVector v_wall_reaction_rate_constant;
    Mutation::Numerics::RealVector v_mass_production_rate;
    std::vector<int> v_species_in_site;

//    double* mp_solution;
    
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const std::vector<Mutation::gsi::CatalysisReaction>& v_reactions;
    
//    WallSteadyStateCoverageSolver* mp_wall_steady_state_coverage_solver;

    double* mp_ropf;
    double* mp_rop;

    GSIStoichiometryManager_FRC v_reactants;
    GSIStoichiometryManager_FRC m_irr_products;

};

    } // namespace gsi
} // namespace Mutation

#endif // CATRATEMANAGER_H
