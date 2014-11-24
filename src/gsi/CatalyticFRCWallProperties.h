#ifndef CATWALLPROPS_H
#define CATWALLPROPS_H

#include "Utilities.h"

namespace Mutation{
    namespace gsi{
      
class CatalysisWallProperties{
    
    CatalysisWallProperties(std::string gsi_wall_properties_file);
    ~CatalysisWallProperties();
    
private:
    int m_nb_active_sites;
    double * m_concentration_active_sites;
    
};
      
    } // namespace gsi
} // namespace Mutation

#endif // CATWALLPROPS_H