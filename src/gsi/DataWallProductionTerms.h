#ifndef DATAWALLPRODUCTIONTERMS_H
#define DATAWALLPRODUCTIONTERMS_H 

#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceDescription.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================

struct DataWallProductionTerms {
    const Mutation::Thermodynamics::Thermodynamics& s_thermo;
    const std::string& s_gsi_mechanism; 
    const Mutation::Utilities::IO::XmlElement& s_node_prod_terms;
    const SurfaceDescription& s_surf_descr;
};

//========================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // DATAWALLPRODUCTIONTERMS_H 
