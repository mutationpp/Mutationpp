/**
 * @file GasSurfaceInteraction.cpp
 *
 * @brief Implements the main class of the GasSurfaceInteraction module.
 */

/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
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
#include "SolidProperties.h"
#include "Surface.h"
#include "SurfaceState.h"

using namespace Eigen;

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
      mp_surf(NULL),
      mp_surf_state(NULL)
{
    if (gsi_input_file == "none"){return;}

    gsi_input_file = databaseFileName(gsi_input_file, "gsi");

    XmlDocument xml_doc(gsi_input_file);
    Mutation::Utilities::IO::XmlElement root_element = xml_doc.root();

    errorWrongTypeofGSIFile(root_element.tag());

    // Improve the error
    root_element.getAttribute("gsi_mechanism", m_gsi_mechanism,
        "For gas-surface interaction a gsi_mechanism should be provided.");

    // Finding the position of the XmlElements
    XmlElement::const_iterator xml_pos_surf_props =
        root_element.findTag("surface_properties");
    XmlElement::const_iterator xml_pos_surf_feats =
        root_element.findTag("surface_features");

    // Creating Surface State class
    mp_surf_state = new SurfaceState(
        m_thermo, *xml_pos_surf_props);

    // Finding the position of the XmlElements
    XmlElement::const_iterator xml_pos_surf_chem =
        root_element.findTag("surface_chemistry");
    XmlElement::const_iterator xml_pos_surf_rad =
        root_element.findTag("surface_radiation");

    // Setting up solid properties
    std::string solid_model;
    if (xml_pos_surf_feats != root_element.end())
        xml_pos_surf_feats->getAttribute("solid_conduction", solid_model);

    XmlElement::const_iterator xml_pos_solid_props;
    if (solid_model == "steady_state") {
        xml_pos_solid_props = root_element.findTag("solid_properties");
        if (xml_pos_solid_props == root_element.end())
            errorSolidPropertiesNotProvided(solid_model);
    } else {
        solid_model = "none";
    }
    DataSolidProperties data_solid_props = { *xml_pos_solid_props };
    mp_surf_state->setSolidProperties(solid_model, data_solid_props);

    // Creating the Surface class
    DataSurface data_surface = {
        m_thermo,
        m_transport,
        m_gsi_mechanism,
        *xml_pos_surf_feats,
        *xml_pos_surf_chem,
        *xml_pos_surf_rad,
        *mp_surf_state };
    mp_surf = Factory<Surface>::create(
        m_gsi_mechanism, data_surface);
}

//==============================================================================

GasSurfaceInteraction::~GasSurfaceInteraction()
{
    if (mp_surf_state != NULL) {delete mp_surf_state;}
    if (mp_surf != NULL) {delete mp_surf;}
}

//==============================================================================

void GasSurfaceInteraction::setSurfaceState(
    const double* const p_mass, const double* const p_energy,
	const int state_variable)
{
    mp_surf_state->setSurfaceState(p_mass, p_energy, state_variable);
}

//==============================================================================

void GasSurfaceInteraction::getSurfaceState(
    double* const p_mass, double* const p_energy,
    const int state_variable)
{
    mp_surf_state->getSurfaceState(p_mass, p_energy, state_variable);
}

//==============================================================================

void GasSurfaceInteraction::surfaceReactionRates(
    double* const p_surf_reac_rates)
{
    VectorXd v_wrk = Map<VectorXd>(p_surf_reac_rates, m_thermo.nSpecies());
    mp_surf->computeSurfaceReactionRates(v_wrk);

    for (int i = 0; i < m_thermo.nSpecies(); i++)
	    p_surf_reac_rates[i] = v_wrk(i);
}

//==============================================================================

void GasSurfaceInteraction::surfaceReactionRatesPerReaction(
    double* const p_surf_rates_per_reac)
{
    Eigen::VectorXd v_surf_rates_per_reac =
        mp_surf->computeSurfaceReactionRatesPerReaction();
	for (int i_r = 0; i_r < mp_surf->nSurfaceReactions(); i_r++){
	    p_surf_rates_per_reac[i_r] = v_surf_rates_per_reac(i_r);
	}
}

//==============================================================================

int GasSurfaceInteraction::nSurfaceReactions()
{
    return mp_surf->nSurfaceReactions();
}

//==============================================================================

void GasSurfaceInteraction::setDiffusionModel(
    const double* const p_mole_frac_edge, const double& dx)
{
    mp_surf->setDiffusionModel(Map<const VectorXd>(
        p_mole_frac_edge, m_thermo.nSpecies()), dx);
}

//==============================================================================

void GasSurfaceInteraction::setGasFourierHeatFluxModel(
    const double* const p_T_edge, const double& dx) {
    mp_surf->setGasFourierHeatFluxModel(Map<const VectorXd>(
        p_T_edge, m_thermo.nEnergyEqns()), dx);
}

//==============================================================================

void GasSurfaceInteraction::setGasRadHeatFlux(
    const double* const m_gas_rad_heat_flux) {
    mp_surf->setGasRadHeatFlux(*m_gas_rad_heat_flux);
}

//==============================================================================

void GasSurfaceInteraction::solveSurfaceBalance()
{
    mp_surf->solveSurfaceBalance();
}

//==============================================================================

void GasSurfaceInteraction::setIterationsSurfaceBalance(const int& iter)
{
    mp_surf->setIterationsSurfaceBalance(iter);
}

//==============================================================================

void GasSurfaceInteraction::getMassBlowingRate(double& mdot){
    mdot = mp_surf->massBlowingRate();
}

//==============================================================================

inline void GasSurfaceInteraction::errorWrongTypeofGSIFile(
    const std::string& gsi_root_tag)
{
    if (gsi_root_tag != "gsi")
    {
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

//==============================================================================
inline void GasSurfaceInteraction::errorSolidPropertiesNotProvided(
    const std::string& error_steady_state)
{
    throw InvalidInputError("GasSurfaceInteraction", error_steady_state)
    << "Solid properties should be provided when steady state assumption "
    << "is assumed for conduction or pyrolysis gases.";
}

    } // namespace GasSurfaceInteraction
} // namespace Mutation
