#include "AutoRegistration.h"
#include "Utilities.h"

#include "DataGSIRateManager.h"
#include "DataGSIReaction.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionsTermsEmpty : public WallProductionTerms {

//======================================================================================

public:
    WallProductionsTermsEmpty( ARGS l_data_wall_prod_terms )
                         : WallProductionTerms(l_data_wall_prod_terms),
                           m_tag("empty") {}

//======================================================================================

    ~WallProductionsTermsEmpty(){}

//======================================================================================

    void productionRate(Eigen::VectorXd& lv_mass_prod_rate){
        lv_mass_prod_rate.setZero();
    }

//======================================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

//======================================================================================
private:
    const std::string m_tag;


};

//======================================================================================

Mutation::Utilities::Config::ObjectProvider<WallProductionsTermsEmpty, WallProductionTerms> wall_productions_terms_empty("empty");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
