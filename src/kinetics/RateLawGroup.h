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

#include "Functors.h"
#include "RateLaws.h"
#include "StateModel.h"

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
     * Returns the temperature used in the last evaluation of the rate
     * coefficients.
     */
    double getT() const { return m_t; }
    
    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk) = 0;
    //{
    //    lnk(p_lnk, Mutation::Numerics::Equals());
    //}
    
    /**
     * Evaluates all of the rates in the group and applies them to the given
     * vector using the given equality operator.
     */
    /*template <typename OP>
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk,
        const OP& op
    ) const = 0;*/

protected:

    /// This is the temperature computed to evaluate the rate law (should be set
    /// in the lnk() function)
    double m_t;
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
    //template <typename OP>
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
        //const OP& op) const
    {
        m_t = TSelectorType().getT(p_state);
        const double lnT  = std::log(m_t);
        const double invT = 1.0 / m_t;
    
        for (int i = 0; i < m_rates.size(); ++i) {
            const std::pair<size_t, RateLawType>& rate = m_rates[i];
            //op(p_lnk[rate.first], rate.second.lnk(lnT));
            p_lnk[rate.first] = rate.second.getLnRate(lnT, invT);
        }
    }

private:

    /// vector of rates to evaluate
    std::vector< std::pair<size_t, RateLawType> > m_rates;
};


/**
 * Manages a collection of RateLawGroup objects such that only one object of any
 * RateLawGroup type is ever created in each collection.
 */
class RateLawGroupCollection
{
public:

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
     * Adds a new rate law to be managed by this collection of rate law groups.
     */
    template <typename GroupType>
    void addRateCoefficient(const size_t rxn, const RateLaw* const p_rate)
    {
        if (m_group_map[&typeid(GroupType)] == NULL)
            m_group_map[&typeid(GroupType)] = new GroupType();
        m_group_map[&typeid(GroupType)]->addRateCoefficient(rxn, p_rate);
        
        std::cout << "Added reaction " << rxn << "("
                  << typeid(GroupType).name() << "), " << "Ngroups = "
                  << m_group_map.size() << std::endl;
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

private:

    /**
     * Small helper class which provides comparison between two std::type_info 
     * pointers.
     */
    struct CompareTypeInfo {
        bool operator ()(const std::type_info* a, const std::type_info* b) const
        {
            return a->before(*b);
        }
    };
    
    typedef std::map<const std::type_info*, RateLawGroup*, CompareTypeInfo>
        GroupMap;
    
    /// Collection of RateLawGroup objects
    GroupMap m_group_map;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // KINETICS_RATE_LAW_GROUP_H

