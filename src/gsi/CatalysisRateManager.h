#ifndef CATRATEMANAGER_H
#define CATRATEMANAGER_H

#include "GSIReaction.h"
#include "GSIStoichiometryManager.h"

/**
 * This class organizes the GSIReactions and computes the necessary matrices and 
 * laws in order to compute in the end the rate law etc. 
 * @todo Fix this description.
 */

namespace Mutation{
    namespace gsi{
      
/**
 * Abstract class for managing the rate of the species for catalytic reactions
 */
class CatalysisRateManager{
public:
    virtual ~CatalysisRateManager(){ }
    virtual CatalysisRateManager* clone() const = 0;

/**
 * This function returns positive production rate if the species are DESTROYED and negative if the species are PRODUCED on the surface.
 * (+) -> Destroyed
 * (-) -> Produced
 */
    virtual void computeRate(const WallState& r_wall_state, double* const p_prod) = 0;
};

class CatalysisGammaRateManager : public CatalysisRateManager{
public:
    /**
     * Default constructor
     */
    CatalysisGammaRateManager( const Mutation::Thermodynamics::Thermodynamics& thermo, const std::vector<CatalysisReaction>& catalytic_reactions );
    
    /**
     * Default destructor
     */
    ~CatalysisGammaRateManager() {
        delete [] mp_gamma;
        delete [] mp_kf;
        delete [] mp_rop;
    }
  
    /**
     * Cloning function
     */
    CatalysisGammaRateManager* clone() const {
        return new CatalysisGammaRateManager(*this);
    }
    
    /**
     * @todo Add description
     */
    void computeRate(const WallState& r_wall_state, double* const p_prod);
    
private:
    /// Number of species in the mixture
    const size_t m_ns;
    
    /// Number of reactions in the catalytic mechanism
    const size_t m_nr;

    double* mp_gamma;
    double* mp_kf;

    double* mp_rop; // DELETE IT

    const std::vector<CatalysisReaction>& m_reactions;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    GSIStoichiometryManager_GammaModel m_reactants;
    GSIStoichiometryManager_GammaModel m_irr_products;
};

class CatalysisFRCRateManager : public CatalysisRateManager{
public:
    /**
     * Default Constructor
     */
    CatalysisFRCRateManager(const Mutation::Thermodynamics::Thermodynamics& thermo,
                            const Mutation::gsi::CatalysisSurfaceProperties* mp_surf_props,
                            const std::vector<Mutation::gsi::CatalysisReaction>& m_catalytic_reactions);
    
    /**
     * Default Destructor
     */
    ~CatalysisFRCRateManager(){ 
        delete [] mp_species_in_site;
        delete [] mp_kf;
        delete [] mp_kb;

        delete [] mp_ropf;
        delete [] mp_rop;
    }

    /**
     * Cloning function
     */
    CatalysisFRCRateManager* clone() const{
        return new CatalysisFRCRateManager(*this);
    }

    void computeRate(const WallState& l_wall_state, double* const p_prod);
    
private:
    void solveSteadyStateCoverage(const WallState& r_wall_state/** @todo is it needed?, double* const p_prod*/);
    
    // Check which of these are actually needed

    const size_t m_ns;
    const size_t m_nr;
    const size_t m_nsites;
    const size_t m_total_species_in_sites;

    double* mp_kf;
    double* mp_kb;
    
    int* mp_species_in_site;
    double* mp_solution;
    
    const std::vector<Mutation::gsi::CatalysisReaction>& m_reactions;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    double* mp_ropf;
    double* mp_ropb;
    double* mp_rop;

    GSIStoichiometryManager_FRC m_reactants;
    GSIStoichiometryManager_FRC m_irr_products;
    // GSIStoichiometryManager m_rev_products;

    // JacobianManagerGSI m_jacobian_manager;

    double* l_prod;
   
};

    } // namespace gsi
} // namespace Mutation

#endif // CATRATEMANAGER_H
