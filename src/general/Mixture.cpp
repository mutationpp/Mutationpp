/**
 * @file Mixture.cpp
 *
 * @brief Mixture class implementation. @see Mutation::Mixture
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "Mixture.h"
#include "StateModel.h"

#include <iostream>
using namespace std;
#include <Eigen/Dense>
using namespace Eigen;

using namespace Mutation::Thermodynamics;

namespace Mutation {

//==============================================================================

Mixture::Mixture(const MixtureOptions& options)
    : Thermodynamics::Thermodynamics(
        options.getSpeciesDescriptor(),
        options.getThermodynamicDatabase(),
        options.getStateModel()),
      Transport(
        *this,
        options.getViscosityAlgorithm(),
        options.getThermalConductivityAlgorithm()),
      Kinetics(
        static_cast<const Thermodynamics&>(*this),
        options.getMechanism()),
      GasSurfaceInteraction(
        *this,
        *this,
        options.getGSIMechanism())
{
    // Add all the compositions given in mixture options to the composition list
    for (int i = 0; i < options.compositions().size(); ++i)
        addComposition(options.compositions()[i]);

    // Set default composition if available
    if (options.hasDefaultComposition())
        setDefaultComposition(m_compositions[options.getDefaultComposition()]);
    
    // Instantiate a new energy transfer model
    state()->initializeTransferModel(*this);
    
}

//==============================================================================

void Mixture::addComposition(const Composition& c, bool make_default)
{
    // Check if all names in the composition are elements
    bool elements = true;
    for (int i = 0; i < c.size(); ++i) {
        if (elementIndex(c[i].name) < 0) {
            elements = false;

            if (speciesIndex(c[i].name) < 0) {
                throw InvalidInputError("composition", c.name())
                    << "Composition has component which is not an element or "
                    << "species belonging to the mixture.";
            }
        }
    }

    // If this composition has species names, then treat all as species and
    // convert to elements
    if (!elements) {
        ArrayXd svals(nSpecies());
        ArrayXd evals(nElements());
        c.getComposition(m_species_indices, svals.data());

        if (c.type() == Composition::MASS)
            convert<Y_TO_X>(svals.data(), svals.data());
        convert<X_TO_XE>(svals.data(), evals.data());

        // Get list of element names
        std::vector<std::string> names;
        for (int i = 0; i < nElements(); ++i)
            names.push_back(elementName(i));

        m_compositions.push_back(
            Composition(names, evals.data(), Composition::MOLE));
    } else
        m_compositions.push_back(c);

    if (make_default)
        setDefaultComposition(m_compositions.back());
}

//==============================================================================

bool Mixture::getComposition(
    const std::string& name, double* const p_vec, Composition::Type type) const
{
    int i = 0;
    for ( ; i < m_compositions.size(); ++i)
        if (m_compositions[i].name() == name) break;

    // Check if there is a composition with the given name
    if (i == m_compositions.size())
        return false;

    // Get the composition
    m_compositions[i].getComposition(m_element_indices, p_vec);

    // Convert if necessary
    if (m_compositions[i].type() != type) {
        if (type == Composition::MOLE)
            convert<YE_TO_XE>(p_vec, p_vec);
        else
            convert<XE_TO_YE>(p_vec, p_vec);
    }

    return true;
}

//==============================================================================
} // namespace Mutation
