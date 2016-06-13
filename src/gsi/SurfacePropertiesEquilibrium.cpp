#include "AutoRegistration.h"
#include "Utilities.h"

#include "Composition.h"
#include "SurfaceProperties.h"

#include <cassert>
#include <iterator>
#include <sstream>

namespace Mutation {
namespace GasSurfaceInteraction {

class SurfacePropertiesEquilibrium : public SurfaceProperties {

public:
    SurfacePropertiesEquilibrium(ARGS l_data_surf_props)
        : SurfaceProperties( l_data_surf_props ),
          idx_sp( l_data_surf_props.s_thermo.nSpecies()),
          idx_ele(l_data_surf_props.s_thermo.nElements()),
          m_condensed_species_descriptor(""),
          m_total_species_descriptor(""),
          m_pyrolysis_blowing( 0.E0 ),
          v_Ykp(l_data_surf_props.s_thermo.nElements()),
          v_Ykc( l_data_surf_props.s_thermo.nElements()),
          m_surface_elemental_constraint() 
    {
        readFromXml(l_data_surf_props.s_node_surf_props, l_data_surf_props.s_thermo);
        mergeSpecies(l_data_surf_props.s_thermo);
    }

    int speciesIndexWall( const std::string& str_sp ) const {return 0;}
    int nSpeciesWall() const {return 0;}

    int nSites() const {return 0;}
    double nTotalSites() const {return 0.0;}

    double fracSite( const int& i_site ) const { return 0.0; }
    int nSpeciesSite( const int& i_site ) const { return 0; }

//======================================================================================

    bool surfaceConstraint() const{
        return m_surface_elemental_constraint;
    }

//======================================================================================

    std::string surfaceSpecies() const{
        return m_total_species_descriptor;
    }

//======================================================================================

    double BprimePyro() const {return m_pyrolysis_blowing;}

//======================================================================================

    void getSurfaceCompositions(
        Eigen::VectorXd& lv_pyrolisysComposition,
        Eigen::VectorXd& lv_charComposition) const
    {
        lv_pyrolisysComposition = v_Ykp;
        lv_charComposition = v_Ykc;
    }

//======================================================================================

    ~SurfacePropertiesEquilibrium(){ }

private:
    std::map<std::string, int> m_element_in;
    bool m_surface_elemental_constraint;
    double m_pyrolysis_blowing;
    std::string m_condensed_species_descriptor;
    std::string m_total_species_descriptor;

    std::vector<Mutation::Thermodynamics::Composition> condensed_elemental_compos;
    int idx_sp;
    int idx_ele;
    Eigen::VectorXd v_Ykp;
    Eigen::VectorXd v_Ykc;

    void mergeSpecies(const Mutation::Thermodynamics::Thermodynamics& surf_thermo){
        std::string recontruct = "";

        for(int i = 0 ; i<idx_sp; ++i){
            recontruct.append(surf_thermo.speciesName(i));
            recontruct.append(" ");
        }

        recontruct.append(m_condensed_species_descriptor);
        Mutation::Utilities::String::trim(recontruct);

        m_total_species_descriptor = recontruct;

    }


    void readFromXml(
        const Mutation::Utilities::IO::XmlElement& surf_props,
        const Mutation::Thermodynamics::Thermodynamics& surf_thermo)
    {

        Mutation::Utilities::IO::XmlElement::const_iterator l_iter_condensed_props;

        for (int i = 0; i < idx_ele; ++i){
            m_element_in[surf_thermo.elementName(i)] = i;
        }

        assert(l_node_surf_props.tag() == "surface_properties");
        // Parse through Environment! Only surface reactions
        Mutation::Utilities::IO::XmlElement::const_iterator l_iter_envs = surf_props.findTag("surface");

        // Checking for surface tags in gsi file
        if (l_iter_envs->tag() == "surface") {
            l_iter_envs->getAttribute("surface_elemental_constraint", m_surface_elemental_constraint, m_surface_elemental_constraint);

            for ( l_iter_condensed_props = l_iter_envs->begin() ;
                  l_iter_condensed_props != l_iter_envs->end() ; l_iter_condensed_props++ ){
                if (l_iter_condensed_props->tag() == "condensed_species"){
                    m_condensed_species_descriptor = Mutation::Utilities::String::trim(l_iter_condensed_props->text());

                }
                else if (l_iter_condensed_props->tag() == "composition"){
                    condensed_elemental_compos.push_back(Mutation::Thermodynamics::Composition(*l_iter_condensed_props));
                    condensed_elemental_compos[condensed_elemental_compos.size()-1].getComposition(m_element_in,&v_Ykc(0));
                }
            }
        }

        l_iter_envs = surf_props.findTag("bulk"); // Fix why when it is commented segmentation fault error from the gsi file
        if (l_iter_envs->tag() == "bulk") {

            l_iter_envs->getAttribute("pyrolysis_mass_blowing_rate", m_pyrolysis_blowing, m_pyrolysis_blowing);

            for ( l_iter_condensed_props = l_iter_envs->begin() ;
                  l_iter_condensed_props != l_iter_envs->end() ; l_iter_condensed_props++ ){
                if (l_iter_condensed_props->tag() == "composition"){
                    condensed_elemental_compos.push_back(Mutation::Thermodynamics::Composition(*l_iter_condensed_props));
                    condensed_elemental_compos[condensed_elemental_compos.size()-1].getComposition(m_element_in, &v_Ykp(0));

                }
            }
        }
        // errorBulkNotFound

    }

}; // class SurfacePropertiesEquilibrium

Mutation::Utilities::Config::ObjectProvider<SurfacePropertiesEquilibrium, SurfaceProperties> surface_properties_equilibrium("equilibrium");

} // namespace GasSurfaceInteraction
} // namespace Mutation 
