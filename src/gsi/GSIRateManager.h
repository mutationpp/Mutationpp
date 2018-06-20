#ifndef GSI_RATE_MANAGER_H
#define GSI_RATE_MANAGER_H

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction;
class SurfaceProperties;
class WallState;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the GSIRateManager class.
 */
struct DataGSIRateManager {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const SurfaceProperties& s_surf_props;
    const WallState& s_wall_state;
    const std::vector<GSIReaction*>& s_reactions;
};

//==============================================================================

/*
 * Abstract class responsible for computing the chemical production rate for
 * Gas-Surface Interaction reactions. Its implementation is different for
 * macroscopic and microscopic approaches.
 */
class GSIRateManager
{
public:
	/**
	 * Required for self registering different rate laws.
	 */
    typedef DataGSIRateManager ARGS;

	/// Returns name of this type.
	static std::string typeName() { return "GSIRateManager"; }

    /**
     * Constructor of the abstract base class, taking as input a structure
     * of type DataGSIRateManager.
     */
    GSIRateManager(ARGS args)
       : m_thermo(args.s_thermo ),
         m_surf_props(args.s_surf_props ),
         m_wall_state(args.s_wall_state ),
         v_reactions(args.s_reactions ) { }

    /**
     * Destructor
     */
    virtual ~GSIRateManager(){ };

    /**
     * This purely virtual function returns the surface production rates
     * according to the selected gas surface interaction model.
     *
     * @return Vector of the gas phase species production rates in kg/m^2-s
     */
    virtual Eigen::VectorXd computeRate() = 0;

protected:

    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceProperties& m_surf_props;
    const WallState& m_wall_state;
    const std::vector<GSIReaction*>& v_reactions;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // GSI_RATE_MANAGER_H
