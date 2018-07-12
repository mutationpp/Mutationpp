/**
 * @file GasSurfaceInteraction.cpp
 *
 * @brief Implements the main class of the GasSurfaceInteraction module.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


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

//==============================================================================

GasSurfaceInteraction::GasSurfaceInteraction(
    Mutation::Thermodynamics::Thermodynamics& thermo,
	Mutation::Transport::Transport& transport,
	std::string gsi_input_file)
    : m_thermo(thermo),
      m_transport(transport),
      mp_surf_solver(NULL),
      mp_surf_props(NULL),
      mp_wall_state(NULL)
{
    if (gsi_input_file == "none"){return;}

    gsi_input_file = databaseFileName(gsi_input_file, "gsi");

    XmlDocument xml_doc(gsi_input_file);
    Mutation::Utilities::IO::XmlElement root_element = xml_doc.root();

    errorWrongTypeofGSIFile(root_element.tag());

    root_element.getAttribute("gsi_mechanism", m_gsi_mechanism, "none");

    // Finding the position of the XmlElements
    Mutation::Utilities::IO::XmlElement::const_iterator xml_pos_surf_props =
        root_element.findTag("surface_properties");
    Mutation::Utilities::IO::XmlElement::const_iterator xml_pos_diff_model =
        root_element.findTag("diffusion_model");
    Mutation::Utilities::IO::XmlElement::const_iterator xml_pos_prod_terms =
        root_element.findTag("production_terms");

    // Creating Surface Properties class
    DataSurfaceProperties data_surface_properties =
        {m_thermo, *xml_pos_surf_props};
    mp_surf_props = Factory<SurfaceProperties>::create(
        m_gsi_mechanism, data_surface_properties);

    // Creating Wall State class
    mp_wall_state = new WallState(m_thermo, *mp_surf_props);

    // Creating the SurfaceBalanceSolver class
    DataSurfaceBalanceSolver data_surface_balance_solver =
        {m_thermo, m_transport, m_gsi_mechanism, *xml_pos_diff_model,
         *xml_pos_prod_terms, *mp_surf_props, *mp_wall_state };
    mp_surf_solver = Factory<SurfaceBalanceSolver>::create(
        m_gsi_mechanism, data_surface_balance_solver );

}

//==============================================================================

GasSurfaceInteraction::~GasSurfaceInteraction()
{
    if (mp_surf_props != NULL) {delete mp_surf_props;}
    if (mp_wall_state != NULL) {delete mp_wall_state;}
    if (mp_surf_solver != NULL) {delete mp_surf_solver;}
}

//==============================================================================

void GasSurfaceInteraction::setWallState(
    const double* const p_mass, const double* const p_energy,
	const int state_variable)
{
    mp_wall_state->setWallState(p_mass, p_energy, state_variable);
}

//==============================================================================

void GasSurfaceInteraction::getWallState(
    double* const p_mass, double* const p_energy,
    const int state_variable)
{
    mp_wall_state->getWallState(p_mass, p_energy, state_variable);
}

//==============================================================================

void GasSurfaceInteraction::surfaceProductionRates(
    double* const p_wall_prod_rates)
{
    Eigen::VectorXd v_wall_rates = mp_surf_solver->computeGSIProductionRates();
	for (int i_sp = 0; i_sp < m_thermo.nSpecies(); i_sp++){
	    p_wall_prod_rates[i_sp] = v_wall_rates(i_sp);
	}
}

//==============================================================================

void GasSurfaceInteraction::setDiffusionModel(
    const double* const p_mole_frac_edge, const double& dx)
{
    mp_surf_solver->setDiffusionModel(Eigen::Map<const Eigen::VectorXd>(
        p_mole_frac_edge, m_thermo.nSpecies()), dx);
}

//==============================================================================

void GasSurfaceInteraction::setConductiveHeatFluxModel(
    const double* const p_T_edge, const double& dx_T)
{
    throw NotImplementedError(
        "GasSurfaceInteraction::setConductiveHeatFluxModel");
}

//==============================================================================

void GasSurfaceInteraction::solveSurfaceBalance()
{
    mp_surf_solver->solveSurfaceBalance();
}

//==============================================================================

void GasSurfaceInteraction::getMassBlowingRate(double& mdot){
    mdot = mp_surf_solver->massBlowingRate();
}

//==============================================================================

void GasSurfaceInteraction::getBprimeCharSpecies(
		std::vector<std::string>& v_species_char_names)
{
    mp_surf_solver->getBprimeCondensedSpecies(v_species_char_names);
}

//==============================================================================

void GasSurfaceInteraction::getBprimeSolution(
    double& bprime_char, std::vector<double>& v_species_char_mass_frac)
{
    mp_surf_solver->getBprimeParameters(bprime_char, v_species_char_mass_frac);
}

//==============================================================================

inline void GasSurfaceInteraction::errorWrongTypeofGSIFile(
    const std::string& gsi_root_tag)
{
    if (gsi_root_tag != "gsi"){
        throw InvalidInputError("GasSurfaceInteraction", gsi_root_tag)
        << "Root element in Gas Surface Interaction input file "
        << gsi_root_tag << " is not of 'gsi' type!";
    }
}

//==============================================================================

inline void GasSurfaceInteraction::errorInvalidGSIFileProperties(
    const std::string& gsi_option)
{
    throw InvalidInputError("GasSurfaceInteraction", gsi_option)
    << gsi_option << " is not a valid gas surface interaction file option!";
}

    } // namespace GasSurfaceInteraction 
} // namespace Mutation
