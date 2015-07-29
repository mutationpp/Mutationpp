#ifndef GSIREACTION_H
#define GSIREACTION_H

#include <string>
#include <vector>

#include "DataGSIReaction.h"
#include "GSIRateLaw.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction {

public:
    typedef const DataGSIReaction& ARGS;

    GSIReaction( ARGS l_data_gsi_reaction )
               : mp_rate_law( NULL ){ }
    virtual ~GSIReaction(){ }

    GSIRateLaw* getRateLaw() const { return mp_rate_law; }

    const std::vector<int>& getReactants() const { return m_reactants; }
    const std::vector<int>& getProducts() const { return m_products; }

protected:
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;

    GSIRateLaw* mp_rate_law;

    inline const char* errorNoFormulainReaction() const { return "No formula specied with reaction!"; }

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSIREACTION_H
