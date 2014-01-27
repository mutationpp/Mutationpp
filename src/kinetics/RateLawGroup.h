/**
 * @file RateLawGroup.h
 *
 * Defines the abstract class RateLawGroup which is used by the RateManager
 * class to combine rate laws which are evaluated at the same temperature in
 * order to compute the rate coefficient for these reactions efficiently.
 *
 * @see class RateLawGroup
 * @see class RateLawGroup1T
 * @see class RateLawGroupCollection
 *
 * @date Oct. 8, 2013
 * @autor J.B. Scoggins (scoggins@vki.ac.be)
 */

#ifndef KINETICS_RATE_LAW_GROUP_H
#define KINETICS_RATE_LAW_GROUP_H

#include <map>
#include <typeinfo>
#include <vector>

#include "RateLaws.h"
#include "Reaction.h"
//#include "StateModel.h"
#include "StoichiometryManager.h"

class StateModel;

namespace Mutation {
    namespace Kinetics {

/**
 * Abstract base class which defines the interface for all RateLawGroup objects
 * which evaluate like rate raws to increase efficiency.
 */
class RateLawGroup
{
public:

    /**
     * Destructor.
     */
    virtual ~RateLawGroup() { };

    /**
     * Adds a new rate to evaluate with this group.
     */
    virtual void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate) = 0;
    
    /**
     * Adds a reaction which uses the backward temperature represented by this
     * group.
     */
    void addReaction(const size_t rxn, const Reaction& reaction) {
        m_reacs.addReaction(rxn, reaction.reactants());
        m_prods.addReaction(rxn, reaction.products());
    }
    
    /**
     * Returns the temperature used in the last evaluation of the rate
     * coefficients.
     */
    double getT() const { return m_t; }
    
    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk) = 0;
        
    /**
     * Computes \Delta G / RT for this rate law group and subtracts these values
     * for each of the reactions in this group.
     */
    void subtractLnKeq(size_t ns, double* const p_g, double* const p_r) const
    {        
        // Compute G_i/RT - ln(Patm/RT)
        const double val = std::log(ONEATM / (RU * m_t));
        for (int i = 0; i < ns; ++i)
            p_g[i] -= val;
        
        // Now subtract \Delta[G_i/RT - ln(Patm/RT)]_j
        m_reacs.decrReactions(p_g, p_r);
        m_prods.incrReactions(p_g, p_r);
        
    }
    

protected:

    /// This is the temperature computed to evaluate the rate law (should be set
    /// in the lnk() function)
    double m_t;
    
    /// Stores the reactants for reactions that will use this rate law for the
    /// reverse direction
    StoichiometryManager m_reacs;
    
    /// Stores the products for reactions that will use this rate law for the
    /// reverse direction
    StoichiometryManager m_prods;
};


/**
 * Groups reaction rate laws based on a single temperature that are the same
 * kind and evaluated at the same temperature together so that they may be 
 * evaluated efficiently.
 */
template <typename RateLawType, typename TSelectorType>
class RateLawGroup1T : public RateLawGroup
{
public:

    /**
     * Adds a new rate to evaluate with this group.
     */
    virtual void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate)
    {
        m_rates.push_back(
            std::make_pair(rxn, dynamic_cast<const RateLawType&>(*p_rate))
        );
    }

    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
    {
        m_t = TSelectorType().getT(p_state);
        const double lnT  = std::log(m_t);
        const double invT = 1.0 / m_t;

        for (int i = 0; i < m_rates.size(); ++i) {
            const std::pair<size_t, RateLawType>& rate = m_rates[i];
            p_lnk[rate.first] = rate.second.getLnRate(lnT, invT);
        }
    }

private:

    /// vector of rates to evaluate
    std::vector< std::pair<size_t, RateLawType> > m_rates;
};


/**
 * Small helper class which provides comparison between two std::type_info
 * pointers.
 */
struct CompareTypeInfo {
    bool operator ()(const std::type_info* a, const std::type_info* b) const {
        return a->before(*b);
    }
};


/**
 * Manages a collection of RateLawGroup objects such that only one object of any
 * RateLawGroup type is ever created in each collection.
 */
class RateLawGroupCollection
{
public:

    typedef std::map<const std::type_info*, RateLawGroup*, CompareTypeInfo>
        GroupMap;

    /**
     * Destructor.
     */
    ~RateLawGroupCollection()
    {
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter) {
            delete iter->second;
            iter->second = NULL;
        }
    }
    
    /**
     * Returns the number of different rate law groups in this collection.
     */
    size_t nGroups() const { return m_group_map.size(); }
    
    /**
     * Returns the GroupMap managed by this RateLawGroupCollection object.
     */
    const GroupMap& groups() const { return m_group_map; }

    /**
     * Adds a new rate law to be managed by this collection of rate law groups.
     */
    template <typename GroupType>
    void addRateCoefficient(const size_t rxn, const RateLaw* const p_rate)
    {
        if (m_group_map[&typeid(GroupType)] == NULL)
            m_group_map[&typeid(GroupType)] = new GroupType();
        m_group_map[&typeid(GroupType)]->addRateCoefficient(rxn, p_rate);
        
        //std::cout << "Added reaction " << rxn << "("
        //          << typeid(GroupType).name() << "), " << "Ngroups = "
        //          << m_group_map.size() << std::endl;
    }
    
    /**
     * Adds a reaction to the manager which allows for the calculation of the 
     * \Delta G / RT term for reverse rate coefficients.
     */
    template <typename GroupType>
    void addReaction(const size_t rxn, const Reaction& reaction)
    {
        if (m_group_map[&typeid(GroupType)] == NULL)
            m_group_map[&typeid(GroupType)] = new GroupType();
        m_group_map[&typeid(GroupType)]->addReaction(rxn, reaction);
    }

    /**
     * Computes the rate coefficients in this collection and stores them in
     * the vector at the index corresponding to their respective reaction.
     */
    void logOfRateCoefficients(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
    {
        // Compute the forward rate constants
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter)
            iter->second->lnk(p_state, p_lnk);
    }
    
    /**
     * Subtracts ln(keq) from the provided rate coefficients.
     */
    void subtractLnKeq(
        const Thermodynamics::Thermodynamics& thermo, double* const p_g,
        double* const p_lnk)
    {
        const size_t ns = thermo.nSpecies();
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter) {
            const RateLawGroup* p_group = iter->second;
            thermo.speciesSTGOverRT(p_group->getT(), p_g);
            p_group->subtractLnKeq(ns, p_g, p_lnk);
        }
    }

private:
    
    /// Collection of RateLawGroup objects
    GroupMap m_group_map;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // KINETICS_RATE_LAW_GROUP_H

