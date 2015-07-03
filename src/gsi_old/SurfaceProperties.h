#ifndef SURFPROPS_H
#define SURFPROPS_H

#include "Utilities.h"
#include "Thermodynamics.h"

namespace Mutation{
    namespace gsi{
      
     class Site{
public:
    Site(const Mutation::Utilities::IO::XmlElement& node, const Mutation::Thermodynamics::Thermodynamics& m_thermo, int& index_species);
    ~Site(){ }
    
    Site( const Site& site )
        : m_number_site_dens(site.m_number_site_dens),
          m_distance(site.m_distance),
          m_label(site.m_label),
          m_gas_species(site.m_gas_species),
          m_wall_species(site.m_wall_species),
          m_str_wall_species(site.m_str_wall_species),
          m_thermo(site.m_thermo) // To remove
    { }
    
    int SitespeciesIndex(const std::string& str_species) const;
    int nSpeciesSite() const { return (m_wall_species.size() - 1);}
    double siteDensity() const { return m_number_site_dens; }

private:
    void parseSiteSpecies(std::string& str, int& index_species);
    
private:
    double m_number_site_dens;
    double m_distance;
    std::string m_label;
    
    std::vector<int> m_gas_species; 
    std::vector<int> m_wall_species; 
    std::vector<std::string> m_str_wall_species;
    
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
};      

class CatalysisSurfaceProperties{
public:
    CatalysisSurfaceProperties(std::string gsi_wall_properties_file, const Mutation::Thermodynamics::Thermodynamics& m_thermo);
    ~CatalysisSurfaceProperties(){
         for( size_t it = 0; it < mp_sites.size() ; ++it){ // Use iterators instead? (const_iter it = mp_sites.begin(); it != mp_site.end(); ++it)
             delete mp_sites[it];
         }
         mp_sites.clear();
    }

    inline int nSites() const { return mp_sites.size(); }
    int nTotalSpeciesinSites() const;
    int nSpeciesinSite(int i_site) const;
    int GSIspeciesIndex(const std::string& str_species) const;
    double getTotalnumberSites() const { return m_total_number_sites; }
    double getnumberinSite(int i_site) const{ return (mp_sites[i_site]->siteDensity() * m_total_number_sites); }

private: 
    void addSite(Site* const site);

    int m_number_envs;
    double m_total_number_sites;
    
    int n_index_sites;
    std::vector<Site*> mp_sites; //Who will destroy them? Access properties? Or interface.

};

//=========================================================================
//=========================================================================
//=========================================================================

class WallState {
public:
    WallState( const Thermodynamics::Thermodynamics& thermo );
    ~WallState(){ }

    void setWallState( const double* const p_rhoi, const double* const p_Twall );
    double* getWallRhoi() const { return p_rhoi_wall; }
    double* getWallNumberDensities() const { return p_number_density_wall; }
    double* getWallTemperature() const { return p_Twall; }
    virtual double* getWallSiteDensities() const { return NULL; } // This is as ugly as it gets.

protected:

    const size_t m_ns;

    Mutation::Numerics::RealVector v_rhoi_wall;
    Mutation::Numerics::RealVector v_number_density_wall;
    double m_Twall;

    double * const p_Twall;
    double * const p_rhoi_wall;
    double * const p_number_density_wall;

    const Thermodynamics::Thermodynamics& m_thermo;

    static const size_t position_of_translational_temperature = 0;

};

//=========================================================================
//=========================================================================
//=========================================================================

class WallStateFRC : public WallState {
public:

    WallStateFRC( const CatalysisSurfaceProperties& surf_props, const Thermodynamics::Thermodynamics& thermo ); 
    ~WallStateFRC(){ }

    double* getWallSiteDensities() const { return p_site_density; }

private:

//    const size_t _nWs; // WHAT IS THIS? WHERE DO I NEED IT?
    const size_t m_number_sites;

    Mutation::Numerics::RealVector v_site_density;
    double * const p_site_density;

    const CatalysisSurfaceProperties& m_surf_props;

};

    } // namespace gsi
} // namespace Mutation

#endif // SURFPROPS_H
