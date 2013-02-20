#include <iostream>
#include <typeinfo>
#include <cmath>

#include "RateManager.h"

using namespace Mutation::Numerics;

namespace Mutation {
    namespace Kinetics {

//==============================================================================

void RateManager::addRateCoefficient(
    const size_t rxn, const RateLaw *const p_rate)
{   
    // Arrhenius
    if (typeid(*p_rate) == typeid(Arrhenius)) {
        m_arrhenius.push_back(
            std::make_pair(rxn, *static_cast<const Arrhenius *const>(p_rate)));
    } else {
        std::cerr << "Rate law " << typeid(*p_rate).name() 
             << " not implemented in RateCoefficientManager!" << std::endl;
        exit(1);
    }
}

//==============================================================================

void RateManager::forwardRateCoefficients(const double T, RealVector& kf) const
{
    lnForwardRateCoefficients(T, kf);
    kf.exp();
}

void RateManager::lnForwardRateCoefficients(
    const double T, RealVector& lnkf) const
{
    std::vector<std::pair<size_t, Arrhenius> >::const_iterator iter;
    
    const double lnT  = log(T);
    const double invT = 1.0 / T;
    
    // Compute all Arrhenius rate laws
    for (iter = m_arrhenius.begin(); iter != m_arrhenius.end(); ++iter)
        lnkf(iter->first) = iter->second.getLnRate(lnT, invT); 
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation

