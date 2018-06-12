#ifndef WALL_PRODUCTION_TERMS_H
#define WALL_PRODUCTION_TERMS_H

#include <eigen3/Eigen/Dense>

namespace Mutation {
    namespace GasSurfaceInteraction {

class Thermodynamics;
class Transport;
class xmlElement;

class GSIReaction;
class GSIRateManager;
class WallState;
class SurfaceProperties;
class WallProductionTerms;

/**
 * Structure which stores the necessary inputs for the
 * WallProductionTerms class.
 */
struct DataWallProductionTerms {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism;
    const Mutation::Utilities::IO::XmlElement& s_node_prod_terms;
    const SurfaceProperties& s_surf_props;
    const WallState& s_wall_state;
    std::vector<WallProductionTerms*>* sp_surf_prod;
    const double* const sp_pres;
};

//==============================================================================

class WallProductionTerms
{
public:
    typedef const DataWallProductionTerms& ARGS;

    WallProductionTerms(ARGS args){}

	/// Returns name of this type.
	static std::string typeName() { return "WallProductionTerms"; }

    virtual ~WallProductionTerms(){}

    virtual void productionRate(Eigen::VectorXd& lv_mass_prod_rate) = 0;

    virtual const std::string& getWallProductionTermTag() const = 0;

}; // class WallProductionTerms

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // WALL_PRODUCTION_TERMS_H
