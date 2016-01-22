#include "AutoRegistration.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

#include <cassert>
#include <iterator>
#include <sstream>

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfacePropertiesFRC : public SurfaceProperties {

public:
SurfacePropertiesFRC( ARGS l_data_surf_props )
                    : SurfaceProperties( l_data_surf_props ),
                      idx_sp( l_data_surf_props.s_thermo.nSpecies() ){

    assert(l_data_surf_props.s_node_surf_props.tag() == "surface_properties" );    

    // Parse through Environment! Only surface reactions 
    Mutation::Utilities::IO::XmlElement::const_iterator l_iter_envs = l_data_surf_props.s_node_surf_props.findTag("surface");
    l_iter_envs->getAttribute("total_number_of_sites", m_tot_site_dens, "No total number of sites for surface has been provided.");
   
    // Parse through Sites
    for (  Mutation::Utilities::IO::XmlElement::const_iterator l_iter_sites = l_iter_envs->begin() ;
           l_iter_sites != l_iter_envs->end() ; l_iter_sites++ ){
        addSite( new Site( *l_iter_sites, l_data_surf_props.s_thermo, idx_sp ) );
        idx_sp++;
    }
    n_sites = vp_sites.size();

}

//=========================================================================
public:
~SurfacePropertiesFRC(){ for( size_t it = 0; it < vp_sites.size() ; ++it) delete vp_sites[it]; }

//=========================================================================

public:
int speciesIndexWall( const std::string& str_sp ) const { // rewrite! // DEBUG!

    int id_sp_wall;
    for ( int i_site = 0; i_site < n_sites; i_site++ ){
        id_sp_wall = vp_sites[i_site]->speciesIndexSite( str_sp );
        if ( id_sp_wall != -1 ) return id_sp_wall;
    }
    return -1;

}

int nSpeciesWall() const {
    int n_species_wall = 0;
    for ( int i_site = 0; i_site < n_sites; i_site++ ){
        n_species_wall += nSpeciesSite( i_site );
    }
    return n_species_wall;
}

int nSites() const { return vp_sites.size(); }
double nTotalSites() const { return m_tot_site_dens; }

double fracSite( const int& i_site ) const { return vp_sites[i_site]->fracSite(); }
int nSpeciesSite( const int& i_site ) const { return vp_sites[i_site]->nSpeciesSite(); }

//===========================================================================================

private:

class Site{
public:
    Site( const Mutation::Utilities::IO::XmlElement& l_node_site, 
          const Mutation::Thermodynamics::Thermodynamics& l_thermo, 
          int& idx_sp ){

        assert( l_node_site.tag() == "site" );

        l_node_site.getAttribute("fraction", m_frac_site_over_surf, 
                                 "The fraction of the surface that is covered with this category of sites should be provided!");
        l_node_site.getAttribute("distance", m_dist, 
                                 "An average distance between these kind of sites should be provided!");
        l_node_site.getAttribute("label", s_label, s_label);

        v_gas_sp.push_back(-1);
        v_wall_sp.push_back(idx_sp); // idx_sp is the global index for all gas and wall species!
        v_str_wall_sp.push_back(s_label);

        std::string l_str_sp = "ErrorinSurfacePropertiesFRC";
        l_node_site.getAttribute("species", l_str_sp, l_str_sp);

        parseSiteSpecies( l_str_sp, l_thermo, idx_sp );

        m_n_species_site = v_wall_sp.size();

    }

    ~Site(){ }

    inline int speciesIndexSite( const std::string& l_str_sp ) const {
        int i_site_sp = 0;
        while( i_site_sp < m_n_species_site ){
            if( l_str_sp == v_str_wall_sp[i_site_sp] ) return v_wall_sp[i_site_sp]; 
            i_site_sp++;
        }
        return -1;
    }

    inline int nSpeciesSite() const { return m_n_species_site; } 

    inline double fracSite() const { return m_frac_site_over_surf; }

private: 
    void parseSiteSpecies( const std::string& s_sp_in_site, const Mutation::Thermodynamics::Thermodynamics& l_thermo , int& idx_sp ){

        std::istringstream iss(s_sp_in_site);
        std::vector<std::string> v_sp_in_site;
        std::copy( std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(v_sp_in_site) );

        for (int i_sp = 0; i_sp < v_sp_in_site.size(); ++i_sp){
            int id_sp = l_thermo.speciesIndex( v_sp_in_site[i_sp] );
            
            if (id_sp == -1){
                std::cerr << "Species " << v_sp_in_site[i_sp] << 
                " in the surface sites are not a species of the gas mixture!" <<
                std::endl;
                exit(1);
            }
            
            v_gas_sp.push_back(id_sp);
            v_str_wall_sp.push_back( v_sp_in_site[i_sp] + '-' + s_label );

            idx_sp++;
            v_wall_sp.push_back(idx_sp);
        }

    }

    double m_frac_site_over_surf;
    double m_dist;
    std::string s_label;

    std::vector<int> v_gas_sp; // Contains a map to the gas species this wall species refers to!
    std::vector<int> v_wall_sp; // Contains the id of the wall species!
    std::vector<std::string> v_str_wall_sp; // Contains the name of the wall species!
    int m_n_species_site;

};

//===========================================================================================

void addSite( Site* const site ){ vp_sites.push_back(site); }
std::vector<Site*> vp_sites;
size_t n_sites;

double m_tot_site_dens;

int idx_sp;

    // Association of some reaction with some kind of surface sites
    

};

Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesFRC, SurfaceProperties> surface_properties_frc("frc");

    } // namespace GasSurfaceInteraction
} // namespace Mutation 
