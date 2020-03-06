/**
 * @file SurfaceState.cpp
 *
 * @brief Class which stores the state of the surface.
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


#include <string>

#include "Errors.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SolidProperties.h"
#include "SurfaceProperties.h"
#include "SurfaceState.h"

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

//==============================================================================

SurfaceState::SurfaceState(
    const Mutation::Thermodynamics::Thermodynamics& thermo,
    const Mutation::Utilities::IO::XmlElement& xml_surf_props)
    : m_thermo(thermo),
      m_ns(thermo.nSpecies()),
      m_nT(thermo.nEnergyEqns()),
      mv_rhoi(m_ns),
      mv_T(m_nT),
      m_is_surface_state_set(false),
      mp_surf_props(NULL),
      mp_solid_props(NULL)
{
    DataSurfaceProperties data_surf_props = {
        m_thermo,
        xml_surf_props };
    std::string name_surf_props;

    xml_surf_props.getAttribute(
        "type", name_surf_props,
        "Error in surface_properties for the gsi imput file. "
        "Surface properties should provided. If no properties are needed "
        "they should set to 'none'.");

    mp_surf_props = Factory<SurfaceProperties>::create(
        name_surf_props, data_surf_props);
}

//==============================================================================

SurfaceState::~SurfaceState(){
    if (mp_surf_props != NULL) {delete mp_surf_props;}
    if (mp_solid_props != NULL) {delete mp_solid_props;}
}

//==============================================================================

void SurfaceState::setSurfaceState(
    const double* const p_mass,
    const double* const p_energy,
    const int state_var)
{
  	switch(state_var){
    case 1:
        setSurfaceRhoi(p_mass);
        setSurfaceT(p_energy);
        break;
    default:
        throw InvalidInputError("variable set", state_var)
        << "This variable-set is not implemented in setSurfaceState"
        << ". Possible variable-sets are:\n"
        << "  1: (species densities, temperature)\n";
    }
    m_is_surface_state_set = true;
}

//==============================================================================

void SurfaceState::getSurfaceState(
    double* const p_rhoi,
    double* const p_rhoie,
    const int state_var) const
{
  	switch(state_var){
    case 1:
        for (int i_sp = 0; i_sp < m_ns; ++i_sp){
            p_rhoi[i_sp] = getSurfaceRhoi()(i_sp);
        }
        for (int i_T = 0; i_T < m_nT ; ++i_T) {
            p_rhoie[i_T] = getSurfaceT()(i_T);
        }
        break;
    default:
        throw InvalidInputError("variable get", state_var)
        << "This variable-get is not implemented in getSurfaceState"
        << ". Possible variable-sets are:\n"
        << "  1: (species densities, temperature)\n";
  	}
}

//==============================================================================

void SurfaceState::setSurfaceRhoi(const double* const p_rhoi){
	mv_rhoi = Eigen::Map<const Eigen::VectorXd>(p_rhoi, m_ns);
}

//==============================================================================

void SurfaceState::setSurfaceT(const double* const p_T){
    mv_T = Eigen::Map<const Eigen::VectorXd>(p_T, m_nT);
}

//==============================================================================

void SurfaceState::setSolidProperties(
    const std::string& solid_model,
    const DataSolidProperties& data_solid_props) {
        mp_solid_props = Factory<SolidProperties>::create(
            solid_model, data_solid_props);
}

    } // GasSurfaceInteraction
} // Mutation
