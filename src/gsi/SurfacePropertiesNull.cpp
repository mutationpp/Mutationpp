#include "AutoRegistration.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfacePropertiesNull : public SurfaceProperties {

public:
    SurfacePropertiesNull( ARGS l_data_surf_props ) : SurfaceProperties( l_data_surf_props ) { }
    ~SurfacePropertiesNull(){ }

    int speciesIndexWall( const std::string& str_sp ) const { return -1; }
    int nSpeciesWall() const { return 0; } 
//    void initializeWallComposition( Eigen::VectorXd& v_wall_comp ) const {}

// private:
//     inline void errorEmptySurfaceProperties(){
//         std::cout << "No Surface Properties have been defined!" << std::endl;
//         exit(1);
//     }

};

Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesNull, SurfaceProperties> surface_properties_null("none");


    } // namespace GasSurfaceInteraction
} // namespace Mutation 
