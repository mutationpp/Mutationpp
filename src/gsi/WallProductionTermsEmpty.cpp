#include "AutoRegistration.h"
#include "Transport.h"
#include "Utilities.h"

#include "WallProductionTerms.h"

using namespace Mutation::Utilities;

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionsTermsEmpty : public WallProductionTerms
{
public:
    WallProductionsTermsEmpty( ARGS l_data_wall_prod_terms )
                         : WallProductionTerms(l_data_wall_prod_terms),
                           m_tag("empty") {}

//==============================================================================

    ~WallProductionsTermsEmpty(){}

//==============================================================================

    void productionRate(Eigen::VectorXd& v_mass_prod_rate){
        v_mass_prod_rate.setZero();
    }

//==============================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

private:
    const std::string m_tag;

};

Config::ObjectProvider<
    WallProductionsTermsEmpty, WallProductionTerms>
    wall_productions_terms_empty("empty");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
