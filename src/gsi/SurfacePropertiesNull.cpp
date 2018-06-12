#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfacePropertiesNull : public SurfaceProperties
{
public:
    SurfacePropertiesNull(ARGS args) : SurfaceProperties(args) { }
    ~SurfacePropertiesNull(){ }

    int speciesIndexWall(const std::string& str_sp) const { return -1; }
    int nSpeciesWall() const { return 0; }

    int nSites() const { return 0; }
    double nTotalSites() const { return 0.0; }

    double fracSite(const int& i_site) const { return 0.0; }
    int nSpeciesSite(const int& i_site) const { return 0; }

}; // class SurfacePropertiesNull

ObjectProvider<
    SurfacePropertiesNull, SurfaceProperties>
    surface_properties_gamma_null("gamma");

    } // namespace GasSurfaceInteraction
} // namespace Mutation 
