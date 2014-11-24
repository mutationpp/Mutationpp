#ifndef CATRATEMANAGER_H
#define CATRATEMANAGER_H

#include "GSIReaction.h"
/**
 * This class organizes the GSIReactions and computes the necessary matrices and 
 * laws in order to compute in the end the rate law etc. 
 * @todo Fix this description.
 */

namespace Mutation{
    namespace gsi{
      
class CatalysisRateManager{
      
public:
    CatalysisRateManager(const std::vector<CatalysisReaction>& v_reaction);
    
      
};

class CatalysisGammaRateManager : CatalysisRateManager{};

    } // namespace gsi
} // namespace Mutation

#endif // CATRATEMANAGER_H