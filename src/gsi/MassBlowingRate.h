#ifndef MASS_BLOWING_RATE_H
#define MASS_BLOWING_RATE_H

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionTerms;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the MassBlowingRate class.
 */
struct DataMassBlowingRate {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    std::vector<WallProductionTerms*>& vs_wall_productions_terms;
};

//==============================================================================

/**
 * Abstract class which returns the mass blowing rate for a heterogeneous
 * reactions. It is the base of a Null object (not used yet) which returns
 * zero blow flux. This is important for the case of purely catalytic or
 * oxidation reactions. In the case where ablative reactions take place,
 * a solid class ablation is available.
 */

class MassBlowingRate
{
public:
	/**
	 * Required for self registering different type of mass blowing rate.
	 */
    typedef const DataMassBlowingRate& ARGS;

	/// Returns name of this type.
	static std::string typeName() { return "MassBlowingRate"; }

    /**
     * Destructor
     */
    virtual ~MassBlowingRate(){ }

    /**
     * Purely virtual function which return the blowing flux in kg/m^2-s.
     */
    virtual double computeBlowingFlux() = 0;

    virtual void getEquilibriumParameters(
        double elementalMassDiff,
        double elementalCharMassFraction,
        double elementalGasMassFraction)
    {
        throw NotImplementedError("MassBlowingRate::getEquilibriumParameters");
    }
};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // MASS_BLOWING_RATE_H
