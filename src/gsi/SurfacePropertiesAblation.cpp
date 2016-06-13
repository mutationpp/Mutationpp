#include "AutoRegistration.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfacePropertiesAblation : public SurfaceProperties
{
public:
    SurfacePropertiesAblation(ARGS args) : SurfaceProperties(args) { }
    ~SurfacePropertiesAblation(){ }

    int speciesIndexWall(const std::string& str_sp) const { return -1; }
    int nSpeciesWall() const { return 0; } 

    int nSites() const { return 0; }
    double nTotalSites() const { return 0.0; }

    double fracSite(const int& i_site) const { return 0.0; }
    int nSpeciesSite(const int& i_site) const { return 0; }

}; // class SurfacePropertiesAblation

Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesAblation, SurfaceProperties> surface_properties_ablation("ablation");

    } // namespace GasSurfaceInteraction
} // namespace Mutation 
