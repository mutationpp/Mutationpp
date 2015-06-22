#include "SurfaceProperties.h"
#include "Utilities.h"
#include <cassert>

using namespace Mutation::Utilities;

namespace Mutation{
    namespace gsi{
      
    /**
     * @todo This should be rewritten for sure. More general for ablation and catalysis case.
     * @todo Check the fraction for the sum is less than 1.
     */

//=========================================================================
    
CatalysisSurfaceProperties::CatalysisSurfaceProperties(std::string gsi_wall_properties_file, const Mutation::Thermodynamics::Thermodynamics& m_thermo)
    : m_total_number_sites(1.E0),
      m_number_envs(0),
      n_index_sites(m_thermo.nSpecies())
{

    gsi_wall_properties_file = getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/gsi/" + gsi_wall_properties_file + ".xml";
    IO::XmlDocument sur(gsi_wall_properties_file);
    IO::XmlElement root = sur.root();

    if(root.tag() != "surface_properties"){
        std::cerr << "Root element in surface properties file " << gsi_wall_properties_file
                      << " is not of 'surface_properties' type!" << std::endl;
        exit(1); 
    }
    /** @todo root.getAttribute("name", name, name); Adding the name of surface properties. No clue if it is actually needed */

    IO::XmlElement::const_iterator iter = root.begin();
    for ( ; iter != root.end(); ++iter) {
        ++m_number_envs;
        if(iter->tag() == "surface"){
            iter->getAttribute("total_number_of_sites", m_total_number_sites, m_total_number_sites);

            IO::XmlElement::const_iterator site = iter->begin();
            for ( ; site !=iter->end(); ++site){
                addSite(new Site(*site, m_thermo, n_index_sites)); 
                ++n_index_sites;
            }
        } else if (iter->tag() == "bulk") {
            std::cerr << "Bulk face not implemented yet!" << std::endl;
            exit(1);
        } else {
            std::cerr << "The tag of the surface properties file should be either surface or bulk" << std::endl;
            exit(1);
        }
    }

}

//=========================================================================

int CatalysisSurfaceProperties::GSIspeciesIndex(const std::string& str_species) const
{
    int i_site = 0;
    int GSIspecies;
    
    while(i_site < nSites()){
        GSIspecies = mp_sites[i_site]->SitespeciesIndex(str_species);
        if(GSIspecies != -1) return GSIspecies;
        i_site++;
    }

    return -1;
}

//=========================================================================

void CatalysisSurfaceProperties::addSite(Site* const site)
{

    mp_sites.push_back(site);

}

//=========================================================================

int CatalysisSurfaceProperties::nTotalSpeciesinSites() const
{

    int total_number_species_in_sites = 0;
    for (int it = 0; it < nSites(); ++it){
        total_number_species_in_sites += mp_sites[it]->nSpeciesSite();
    }
    return total_number_species_in_sites;
  
}

//=========================================================================

int CatalysisSurfaceProperties::nSpeciesinSite ( int i_site ) const
{

    return mp_sites[i_site]->nSpeciesSite();

}

//=========================================================================
//=========================================================================
//=========================================================================

Site::Site( const IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& thermo, int& index_species)
          : m_thermo(thermo)
{

Mutation::Utilities::IO::XmlElement::const_iterator iter = node.begin();

    assert(node.tag() == "site");

        node.getAttribute("fraction", m_number_site_dens, m_number_site_dens);
        node.getAttribute("distance", m_distance, m_distance);
        node.getAttribute("label", m_label, m_label);

        m_gas_species.push_back(-1);
        m_wall_species.push_back(index_species); // This might not be actually needed.
        m_str_wall_species.push_back(m_label);
        
        std::string m_string_species;
        node.getAttribute("species", m_string_species, m_string_species);

        m_string_species = String::trim(m_string_species);
        parseSiteSpecies(m_string_species, index_species);

}

//=========================================================================

int Site::SitespeciesIndex( const std::string& str_species ) const
{

    int i_site_species = 0;

    while(i_site_species < m_str_wall_species.size()){
        if(str_species == m_str_wall_species[i_site_species]) return m_wall_species[i_site_species]; /** @todo Maybe a better comparison */
        i_site_species++;
    }
    
    return -1;
}

//=========================================================================

void Site::parseSiteSpecies( std::string& str, int& index_species ){

    /**
     * @todo What happens if I have two blank spaces between the species? -It gives wrong results 
     */

    size_t c = 0;
    size_t s = 0;
    int index;

    while(c != str.size()){
        if(str[c] == ' ' && s != 0){
            index = m_thermo.speciesIndex(str.substr(c-s,c));
            if (index == -1){
                std::cerr << "Species " << str.substr(c-s,c) << " in the surface sites are not a species of the gas mixture!" << std::endl;
                exit(1);
            }
            m_gas_species.push_back(index);
            ++index_species;
            m_wall_species.push_back(index_species);
            m_str_wall_species.push_back(str.substr(c-s,c) + "-" + m_label);
            s = 0;
        } else s++;
        c++;
    }
        index = m_thermo.speciesIndex(str.substr(c-s,c));
        if (index == -1){
            std::cerr << "Species " << str.substr(c-s,c) << " in the surface sites are not a species of the gas mixture!" << std::endl;
            exit(1);
        }
        m_gas_species.push_back(index);
        m_str_wall_species.push_back(str.substr(c-s,c) + "-" + m_label);
        ++index_species;
        m_wall_species.push_back(index_species);
        
}

//=========================================================================
//=========================================================================
//=========================================================================

/**
 * @todo Possible bug in  line:
 * _sitedens[n_counter] = l_site_dens * ( 1.E0 -  static_cast<double>((n_species_in_site)-1) * inv_species_in_site);
 */

WallState::WallState( const Thermodynamics::Thermodynamics& thermo )
                      : m_ns(thermo.nSpecies()),
                        v_rhoi_wall(thermo.nSpecies()),
                        v_number_density_wall(thermo.nSpecies()),
                        m_thermo(thermo),
                        p_Twall(&m_Twall),
                        p_rhoi_wall( &v_rhoi_wall(0) ),
                        p_number_density_wall( &v_number_density_wall(0) )
{ }

//=========================================================================

void WallState::setWallState(const double* const lp_rhoi, const double* const lp_Twall){

    for (int i_ns = 0; i_ns < m_ns; ++i_ns){
        v_rhoi_wall(i_ns) = lp_rhoi[i_ns];
        v_number_density_wall(i_ns) = v_rhoi_wall(i_ns)/(m_thermo.speciesMw(i_ns)) * Mutation::NA ; 
    }

    m_Twall = lp_Twall[position_of_translational_temperature];

}
//=========================================================================
//=========================================================================
//=========================================================================

/**
 * @todo Possible bug in  line:
 * _sitedens[n_counter] = l_site_dens * ( 1.E0 -  static_cast<double>((n_species_in_site)-1) * inv_species_in_site);
 */

WallStateFRC::WallStateFRC( const CatalysisSurfaceProperties& surf_props, const Thermodynamics::Thermodynamics& thermo )
                      : WallState( thermo ),
                        m_surf_props( surf_props ),
                        m_number_sites( surf_props.nSites() ),
                        v_site_density( surf_props.nTotalSpeciesinSites() + surf_props.nSites() ), 
                        p_site_density( &v_site_density(0) )
{
    /** @todo This is not designed well. */

    int n_species_in_site = 0;
    int n_counter = 0;
    double l_site_dens = 0.E0;
    double inv_species_in_site = 0.E0;

    for (int i_sites = 0; i_sites < m_number_sites; ++i_sites){

        n_species_in_site = surf_props.nSpeciesinSite(i_sites);
        inv_species_in_site  = 1.E0 / static_cast<double>(n_species_in_site);
        l_site_dens =  surf_props.getnumberinSite(i_sites);
        v_site_density(n_counter) = l_site_dens * ( 1.E0 - static_cast<double>((n_species_in_site)-1) * inv_species_in_site);
        ++n_counter;
        for (int i_species_in_site = 1; i_species_in_site < n_species_in_site + 1; ++i_species_in_site){
            v_site_density(n_counter) = l_site_dens * inv_species_in_site;
            ++n_counter;
        }
    }

}

//==========================================================================

    } // namespace gsi
} // namespace Mutation
