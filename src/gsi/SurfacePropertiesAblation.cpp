/**
 * @file SurfacePropertiesAblation.cpp
 *
 * @brief SurfaceProperties class when an ablative surface is modeled
 *        based on a gamma model.
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

#include <iterator>

#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

using namespace std;

using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Class which stores the properties necessary for the modelling of
 * an ablative surface. The surface is supposed to be composed of smaller
 * surfaces, like solid carbon or resin. Each different pseudo surfaces is
 * distinguished by a different label.
 */
class SurfacePropertiesAblation : public SurfaceProperties
{
public:
	SurfacePropertiesAblation(ARGS args)
        : SurfaceProperties(args),
          xml_surf_props(args.s_node_surf_props),
          m_thermo(args.s_thermo),
          n_gas_sp(m_thermo.nSpecies()),
          is_surface_set(false)
    {
        assert(xml_surf_props.tag() == "surface_properties");

        // For all bulk compositions
        for (XmlElement::const_iterator iter_phase =
                xml_surf_props.begin();
            iter_phase != xml_surf_props.end();
            iter_phase++)
        {
            string option = iter_phase->tag();

            if (option.compare("surface") == 0) {
                std::string label;
                iter_phase->getAttribute(
                    "label", label,
                    "Error in surface option for the "
                    "surface properties. A label should be provided for "
                    "this type of surface.");

                std::string species;
                iter_phase->getAttribute(
                    "species", species,
                    "Error in surface option for the "
                    "surface properties. Species should be provided for "
                    "this type of surface.");
                parseAblationSpecies(species, label);

                is_surface_set = true;
            } else {
                throw InvalidInputError("SurfaceProperties", option)
                 << option << "is a wrong input for surface "
                << "properties.";
            }

            n_surf_sp = v_surf_sp.size();
        }

        if (!is_surface_set){
            throw InvalidInputError("SurfaceProperties", xml_surf_props.tag())
            << "In the surface properties at least one type of surface "
            << "should be provided.";
        }

    }

//==============================================================================
    /**
     * Destructor.
     */
    ~SurfacePropertiesAblation(){ }

//==============================================================================

    /**
     * Returns the index of the surface species according to the order they
     * appear in the input file, following the gas phase species.
     */
     int surfaceSpeciesIndex(const std::string& str_sp) const {
        for (int i_sp = 0; i_sp < n_surf_sp; i_sp++) {
            if (v_surf_sp[i_sp] == str_sp)
                return n_gas_sp + i_sp;
        }
        return -1;
    }

//==============================================================================
    /**
     * Returns the gas phase species associated with the surface species.
     */
    int surfaceToGasIndex(const int& i_surf_species) const {
        return v_surf_to_gas_idx[i_surf_species - n_gas_sp];
    }

//==============================================================================
    /**
     * Returns the number of surface species.
     */
    size_t nSurfaceSpecies() const { return n_surf_sp; }

//==============================================================================
private:
    void parseAblationSpecies(
        const std::string& species,
        const std::string& label)
    {
        std::istringstream iss(species);
        std::vector<std::string> v_species;
        std::copy(std::istream_iterator<std::string>(iss),
                  std::istream_iterator<std::string>(),
                  std::back_inserter(v_species));

        for (int i_sp = 0; i_sp < v_species.size(); ++i_sp)
        {
            int id_sp = m_thermo.speciesIndex(v_species[i_sp]);

            if (id_sp == -1) {
                throw InvalidInputError("SurfaceProperties",
                    v_species[i_sp]) << "Surface species " <<
                    v_species[i_sp] << " is not " <<
                    "a species of the gas mixture!";
            }

            v_surf_to_gas_idx.push_back(id_sp);
            v_surf_sp.push_back(v_species[i_sp] + '-' + label);
        }
    }

private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const XmlElement& xml_surf_props;

    const size_t n_gas_sp;
    int n_surf_sp;

    bool is_surface_set;

    std::vector<std::string> v_surf_sp;
    std::vector<int> v_surf_to_gas_idx;
};

ObjectProvider<
    SurfacePropertiesAblation, SurfaceProperties>
    surface_properties_ablation("ablation");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
