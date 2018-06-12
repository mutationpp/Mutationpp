#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "GasSurfaceInteraction.h"
#include "SurfaceBalanceSolver.h"
#include "SurfaceProperties.h"
#include "WallState.h"

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

GasSurfaceInteraction::GasSurfaceInteraction(
    Mutation::Thermodynamics::Thermodynamics& l_thermo,
	Mutation::Transport::Transport& l_transport,
	std::string l_gsi_input_file)
    : m_thermo(l_thermo),
      m_transport(l_transport),
      mp_surf_solver(NULL),
      mp_surf_props(NULL),
      mp_wall_state(NULL)
{
    if (l_gsi_input_file == "none"){return;}

    l_gsi_input_file = databaseFileName(l_gsi_input_file, "gsi");

    XmlDocument l_xml_doc(l_gsi_input_file);
    XmlElement l_root_element = l_xml_doc.root();

    errorWrongTypeofGSIFile(l_root_element.tag());

    l_root_element.getAttribute("gsi_mechanism", m_gsi_mechanism, "none");

    // Finding the position of the XmlElements
    XmlElement::const_iterator xml_pos_surf_props =
        l_root_element.findTag("surface_properties");
    XmlElement::const_iterator xml_pos_diff_model =
        l_root_element.findTag("diffusion_model");
    XmlElement::const_iterator xml_pos_prod_terms =
        l_root_element.findTag("production_terms");

    // Creating Surface Properties class
    // xml_pos_surf_props.tag("none" or whatever) and then m_gsi_mechanism->to this
    DataSurfaceProperties l_data_surface_properties =
        {m_thermo, *xml_pos_surf_props};
    mp_surf_props = Factory<SurfaceProperties>::create(
    		m_gsi_mechanism, l_data_surface_properties);

    // Creating Wall State class
    mp_wall_state = new WallState(m_thermo, *mp_surf_props);

    // Creating the SurfaceBalanceSolver class
    DataSurfaceBalanceSolver l_data_surface_balance_solver =
        {m_thermo, m_transport, m_gsi_mechanism, *xml_pos_diff_model,
         *xml_pos_prod_terms, *mp_surf_props, *mp_wall_state };
    mp_surf_solver = Factory<SurfaceBalanceSolver>::create(
        m_gsi_mechanism, l_data_surface_balance_solver );

}

//======================================================================================

GasSurfaceInteraction::~GasSurfaceInteraction()
{
    if (mp_surf_props != NULL) {delete mp_surf_props;}
    if (mp_wall_state != NULL) {delete mp_wall_state;}
    if (mp_surf_solver != NULL) {delete mp_surf_solver;}
}

//======================================================================================

void GasSurfaceInteraction::setWallState(
    const double* const p_mass, const double* const p_energy,
	const int state_variable)
{
    mp_wall_state->setWallState(p_mass, p_energy, state_variable);
}

//======================================================================================

void GasSurfaceInteraction::getWallState(
    double* const p_mass, double* const p_energy,
    const int state_variable)
{
    mp_wall_state->getWallState(p_mass, p_energy, state_variable);
}

//======================================================================================

void GasSurfaceInteraction::surfaceProductionRates(double* const p_wall_prod_rates)
{
    Eigen::VectorXd v_wall_rates = mp_surf_solver->computeGSIProductionRates();
	for (int i_sp = 0; i_sp < m_thermo.nSpecies(); i_sp++){
	    p_wall_prod_rates[i_sp] = v_wall_rates(i_sp);
	}
}

//======================================================================================

void GasSurfaceInteraction::setDiffusionModel(
    const double* const p_mole_frac_edge, const double& dx)
{
    mp_surf_solver->setDiffusionModel(Eigen::Map<const Eigen::VectorXd>(
        p_mole_frac_edge, m_thermo.nSpecies()), dx);
}

//======================================================================================

void GasSurfaceInteraction::setConductiveHeatFluxModel( // Experimental
    const double* const p_T_edge, const double& dx_T)
{
//    mp_surf_solver->setConductiveHeatFluxModel(p_T_edge, dx);
}

//======================================================================================

void GasSurfaceInteraction::solveSurfaceBalance()
{
    mp_surf_solver->solveSurfaceBalance();
}

//=================================================================

void GasSurfaceInteraction::getMassBlowingRate(double& mdot){
    mdot = mp_surf_solver->massBlowingRate();
}

//=================================================================

void GasSurfaceInteraction::getBprimeCharSpecies(
		std::vector<std::string>& v_species_char_names)
{
    mp_surf_solver->getBprimeCondensedSpecies(v_species_char_names);
}

//=================================================================

void GasSurfaceInteraction::getBprimeSolution(
    double& bprime_char, std::vector<double>& v_species_char_mass_frac)
{
    mp_surf_solver->getBprimeParameters(bprime_char, v_species_char_mass_frac);
}

//======================================================================================

inline void GasSurfaceInteraction::errorWrongTypeofGSIFile(const std::string& gsi_root_tag)
{
    if (gsi_root_tag != "gsi"){
        std::cerr << "Root element in Gas Surface Interaction input file " << gsi_root_tag
        		  << " is not of 'gsi' type!" << std::endl; // @todo FIX ERROR not l_gsi_root_tag. Instead name of file...
        exit(1);
    }
}

//======================================================================================

inline void GasSurfaceInteraction::errorInvalidGSIFileProperties(const std::string& gsi_option)
{
    std::cerr << gsi_option << " is not a valid gas surface interaction file option!"
              << std::endl;
    exit(1);
}

//======================================================================================

    } // namespace GasSurfaceInteraction 
} // namespace Mutation
