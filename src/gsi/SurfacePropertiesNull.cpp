#include "AutoRegistration.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

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

Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesNull, SurfaceProperties> surface_properties_gamma_null("gamma");
Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesNull, SurfaceProperties> surface_properties_gamma_energy_null("gamma_energy");

    } // namespace GasSurfaceInteraction
} // namespace Mutation 
