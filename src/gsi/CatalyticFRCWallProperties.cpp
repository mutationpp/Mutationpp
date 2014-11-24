#include "CatalyticFRCWallProperties.h"
#include "Utilities.h"

using namespace Mutation::Utilities;

/**
* This is a class add description. Reads and stores properties from Catalysis_Nitrogen_Silica_Rec...
*/

namespace Mutation{
    namespace gsi{

CatalysisWallProperties::CatalysisWallProperties(std::string gsi_wall_properties_file):
                                        m_nb_active_sites(0)
{
    gsi_wall_properties_file = 
    getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/gsi/" +
    gsi_wall_properties_file + ".xml";

    IO::XmlDocument doc(gsi_wall_properties_file);        
    IO::XmlElement root = doc.root();

    if (root.tag() != "FRC_data") {
        std::cout << "Root element in gsi_wall_properties_file " << gsi_wall_properties_file
                  << " is not of 'FRC_data' type!";
        exit(1); 
    }
    
    m_concentration_active_sites = new double [1];

}
      
CatalysisWallProperties::~CatalysisWallProperties()
{
delete [] m_concentration_active_sites;
}

      
    } // namespace gsi
} // namespace Mutation