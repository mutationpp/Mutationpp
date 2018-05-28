#include "AutoRegistration.h"
#include "Utilities.h"
#include "Composition.h"

#include "DataGSIRateManager.h"
#include "DataGSIReaction.h"
#include "WallProductionTerms.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

class WallProductionSurfaceChemistry;

class WallProductionPyrolysis : public WallProductionTerms
{
public:
    WallProductionPyrolysis(ARGS args)
        : WallProductionTerms(args),
          m_thermo(args.s_thermo),
          m_ns(m_thermo.nSpecies()),
          m_neqns(m_ns),
          m_wall_state(args.s_wall_state),
          m_tag("pyrolysis"),
          m_mass_char(0.0),
          v_el_comp(m_thermo.nElements()),
          v_equil_comp(m_ns),
          p_Pwall(args.sp_pres)
    {
         Mutation::Utilities::IO::XmlElement::const_iterator iter_xml_elem_data;

         // Getting the data
         iter_xml_elem_data = args.s_node_prod_terms.findTag("properties");

         iter_xml_elem_data->getAttribute("phi", m_phi, "The ratio of the charred and virgin material should be provided.");

         std::map<std::string, int> l_map;
         for(int i_elem = 0; i_elem < m_thermo.nElements(); i_elem++){
             l_map[m_thermo.elementName(i_elem)] = i_elem;
         }

         Mutation::Utilities::IO::XmlElement xml_comp(*args.s_node_prod_terms.findTag("composition"));

         Mutation::Thermodynamics::Composition m_comp(xml_comp);
         m_comp.getComposition(l_map, v_el_comp.data());

         pv_wall_prod = args.sp_surf_prod;
    }

//======================================================================================

    ~WallProductionPyrolysis(){}

//======================================================================================

    void productionRate(Eigen::VectorXd& lv_pyrolysis_mass_source)
    {
        lv_pyrolysis_mass_source.setZero();
        static Eigen::VectorXd v_wrk(m_neqns);

        m_mass_char = 0.0;
        for(int i_term = 0; i_term < pv_wall_prod->size(); ++i_term)
        {
        	v_wrk.setZero();
            if(((*pv_wall_prod)[i_term]->getWallProductionTermTag()).compare("surface_chemistry") == 0){
                (*pv_wall_prod)[i_term]->productionRate(v_wrk);
                m_mass_char += v_wrk.sum();
            }
        }

        const size_t pos_T_trans = 0;
        double Twall = m_wall_state.getWallT()(pos_T_trans);

        m_thermo.equilibriumComposition(Twall, *p_Pwall,
                                        v_el_comp.data(),
                                        v_equil_comp.data());

        // Units of v_equil_comp
        m_thermo.convert<Mutation::Thermodynamics::X_TO_Y>(v_equil_comp.data(), v_equil_comp.data());

        lv_pyrolysis_mass_source = v_equil_comp * m_phi * m_mass_char;

    }

//======================================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

//======================================================================================
private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    const WallState& m_wall_state;
    std::vector<WallProductionTerms*>* pv_wall_prod;

    const size_t m_ns;
    const size_t m_neqns;

    const std::string m_tag;

    double m_phi;
    double m_mass_char;

    const double* const p_Pwall;

    Eigen::VectorXd v_el_comp;
    Eigen::VectorXd v_equil_comp;

//======================================================================================

}; // class SurfaceProductionsRates

//======================================================================================

Mutation::Utilities::Config::ObjectProvider<WallProductionPyrolysis, WallProductionTerms> wall_production_pyrolysis("pyrolysis");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
