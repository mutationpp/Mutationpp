#ifndef KINETICS_RATEMANAGER_H
#define KINETICS_RATEMANAGER_H

#include "RateLawGroup.h"
#include "Reaction.h"

namespace Mutation {
    namespace Kinetics {

class Reaction;

/**
 * Manages the efficient computation of rate coefficients for the evaluation of 
 * reaction rates.
 */
class RateManager
{
public:

    /**
     * Creates a new RateManager which manages the reaction rate coefficients
     * for the given reactions.
     */
    RateManager(size_t ns, const std::vector<Reaction>& reactions);
    
    /**
     * Destructor.
     */
    ~RateManager();
    
    /**
     * Updates the current values of the rate coefficients.
     */
    void update(const Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns a pointer to the forward rate coefficients evaluated at the 
     * forward temperature.
     */
    const double* const lnkf() { return mp_lnkf; }
    
    /**
     * Returns a pointer to the forward rate coefficients evaluated at the 
     * forward temperature.
     */
    const double* const lnkb() { return mp_lnkb; }

private:

    /**
     * Add a new reaction whose rate law should be managed.  This method will
     * select the appropriate temperature to evaluate the given reaction's rate
     * law based on the type of reaction it is.
     */
    void addReaction(const size_t rxn, const Reaction& reactions);

    /**
     * Helper function which enumerates all of the possible reaction types and
     * calls the appropriate addRate() template based on the reaction type.
     */
    template <int NReactionTypes>
    void selectRate(const size_t rxn, const Reaction& reaction);

    /**
     * Small helper function which adds a new rate law to the manager based on
     * which rate law groups the forward and reverse rates fall into.
     */
    template <typename ForwardGroup, typename ReverseGroup>
    void addRate(const size_t rxn, const Reaction& reaction);

private:

    /// Number of species in the mixture
    const size_t m_ns;
    
    /// Number of reactions in the mechanism
    const size_t m_nr;

    /// Collection of different rate law groups (can be forward or reverse)
    RateLawGroupCollection m_rate_groups;
    
    /// Storage for forward rate coefficient at forward temperature
    double* mp_lnkf;
    
    /// Storage for forward rate coefficient at backward temperature
    double* mp_lnkb;
    
    /// Storage for species Gibbs free energies
    double* mp_gibbs;
    
    /// Stores the indices for which the forward and reverse temperature
    /// evaluations are equal
    std::vector<size_t> m_to_copy;
};


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_RATEMANAGER_H
