#include "AutoRegistration.h"
#include "Composition.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "WallProductionTerms.h"
#include "WallState.h"

#include "GSIRateManager.h"
#include "GSIReaction.h"

using namespace Eigen;

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace GasSurfaceInteraction {

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
          mv_el_comp(m_thermo.nElements()),
          mv_equil_comp(m_ns),
          mp_Pwall(args.sp_pres)
    {
         XmlElement::const_iterator
             iter_xml_elem_data;

         // Getting the data
         iter_xml_elem_data = args.s_node_prod_terms.findTag("properties");

         iter_xml_elem_data->getAttribute("phi", m_phi,
             "The ratio of the charred and virgin material phi "
             "should be provided.");

         std::map<std::string, int> map;
         for(int i_elem = 0; i_elem < m_thermo.nElements(); i_elem++) {
             map[m_thermo.elementName(i_elem)] = i_elem;
         }

         XmlElement xml_comp(
             *args.s_node_prod_terms.findTag("composition"));

         Composition m_comp(xml_comp);
         m_comp.getComposition(map, mv_el_comp.data());

         mpv_wall_prod = args.sp_surf_prod;
    }

//==============================================================================

    ~WallProductionPyrolysis(){}

//==============================================================================

    void productionRate(VectorXd& v_pyrolysis_mass_source)
    {
        v_pyrolysis_mass_source.setZero();
        static VectorXd v_wrk(m_neqns);

        // At steady state the mass of the pyrolysis gases are proportional to
        // the mass of the blowing gases. The following part computes the mass
        // produced due to ablation reactions. From all the surface terms, only
        // the ones with the tag "surface_chemistry" are considered.
        m_mass_char = 0.0;
        for(int i_term = 0; i_term < mpv_wall_prod->size(); ++i_term)
        {
        	v_wrk.setZero();
            if(((*mpv_wall_prod)[i_term]->
                getWallProductionTermTag()).compare("surface_chemistry") == 0)
            {
                (*mpv_wall_prod)[i_term]->productionRate(v_wrk);
                m_mass_char += v_wrk.sum();
            }
        }

        const size_t pos_T_trans = 0;
        double Twall = m_wall_state.getWallT()(pos_T_trans);

        m_thermo.equilibriumComposition(Twall, *mp_Pwall,
                                        mv_el_comp.data(),
                                        mv_equil_comp.data());

        // Units of v_equil_comp
        m_thermo.convert<X_TO_Y>(mv_equil_comp.data(), mv_equil_comp.data());

        v_pyrolysis_mass_source = mv_equil_comp * m_phi * m_mass_char;

    }

//==============================================================================

    const std::string& getWallProductionTermTag() const { return m_tag; }

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    const WallState& m_wall_state;
    std::vector<WallProductionTerms*>* mpv_wall_prod;

    const size_t m_ns;
    const size_t m_neqns;

    const std::string m_tag;

    double m_phi;
    double m_mass_char;

    const double* const mp_Pwall;

    VectorXd mv_el_comp;
    VectorXd mv_equil_comp;

};

Config::ObjectProvider<
    WallProductionPyrolysis, WallProductionTerms>
    wall_production_pyrolysis("pyrolysis");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
