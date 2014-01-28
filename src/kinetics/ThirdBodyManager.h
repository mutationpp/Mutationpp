#ifndef KINETICS_THIRDBODYMANAGER_H
#define KINETICS_THIRDBODYMANAGER_H

#include <vector>
#include <utility>

namespace Mutation {
    namespace Kinetics {

/**
 * Small helper class for class ThirdbodyManager which a thirdbody update on a
 * single reaction.
 */
class PartialThirdbodyEffs
{
public:

    PartialThirdbodyEffs(
        const size_t rxn, const std::vector<std::pair<int, double> >& effs)
        : m_rxn(rxn), m_effs(effs)
    { }
    
    inline void multiplyEfficiencies(
        double sum, const double* const p_s, double* const p_r) const
    {
        std::vector<std::pair<int, double> >::const_iterator iter;
        for (iter = m_effs.begin(); iter != m_effs.end(); ++iter)
            sum += p_s[iter->first] * iter->second;
        p_r[m_rxn] *= sum;
    }
    
private:
    
    size_t m_rxn;
    std::vector<std::pair<int, double> > m_effs;
    
}; // class PartialThirdbodyEffs


/**
 * Manages the efficient application of thirdbody terms to reaction rates of 
 * progress.
 */
class ThirdbodyManager
{
public:

    /**
     * Constructor
     */
    ThirdbodyManager(const size_t ns)
        : m_ns(ns)
    { }
    
    /**
     * Adds a new thirdbody reaction to be managed by this manager.
     */
    void addReaction(
        const size_t rxn, const std::vector<std::pair<int, double> > effs)
    {
        std::vector<std::pair<int, double> > partial_effs;
        std::vector<std::pair<int, double> >::const_iterator iter;
        
        for (iter = effs.begin(); iter != effs.end(); ++iter) {
            if (iter->second != 1.0)
                partial_effs.push_back(
                    std::make_pair(iter->first, iter->second - 1.0));
        }
        
        m_effs.push_back(PartialThirdbodyEffs(rxn, partial_effs));
    }

    /**
     * Multiplies the thirdbody reaction rates of progress by their 
     * corresponding thirdbody efficiency sums given the species molar 
     * concentrations vector.
     */
    void multiplyThirdbodies(const double* const p_s, double* const p_r) const
    {
        double sum = 0.0;
        for (int i = 0; i < m_ns; ++i)
            sum += p_s[i];
        
        std::vector<PartialThirdbodyEffs>::const_iterator iter = m_effs.begin();
        for ( ; iter != m_effs.end(); ++iter)
            iter->multiplyEfficiencies(sum, p_s, p_r);
    }

private:

    const size_t m_ns;
    std::vector<PartialThirdbodyEffs> m_effs;
    
}; // class ThirdbodyManager


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_THIRDBODYMANAGER_H
