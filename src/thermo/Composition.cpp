/**
 * @file Composition.cpp
 *
 * @brief Implementation of the Thermodynamics::Composition class.
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

#include <iostream>
#include <cstdlib>

#include "Composition.h"
#include "Errors.h"
#include "StringUtils.h"
#include "XMLite.h"

using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//=============================================================================

Composition::Composition(const char* const list)
    : m_name(""), m_type(MOLE)
{
    std::string error = componentsFromList(list);
    if (error != "")
        throw InvalidInputError("composition", list) << error;
}

//=============================================================================

Composition::Composition(
    const std::string& name, const std::string& list, Composition::Type type)
    : m_name(name), m_type(type)
{
    std::string error = componentsFromList(list);
    if (error != "")
        throw InvalidInputError("composition", list)("name", name) << error;
}

//=============================================================================

Composition::Composition(const IO::XmlElement& element)
{
    // Get the composition name
    element.getAttribute("name", m_name, "A composition must have a name!");

    // Get the type of composition (mass or mole)
    std::string type_string("mole");
    element.getAttribute("type", type_string, type_string);
    type_string = String::toLowerCase(type_string);

    if (type_string == "mole")
        m_type = MOLE;
    else if (type_string == "mass")
        m_type = MASS;
    else
        element.parseError("Composition type can be either 'mass' or 'mole'!");

    std::string error = componentsFromList(element.text());
    if (error != "")
        element.parseError(error);
}

//=============================================================================

Composition::Composition(
    const std::vector<std::string>& names, const double* const vals, Type type)
    : m_type(type)
{
    for (int i = 0; i < names.size(); ++i)
        if (vals[i] != 0.0)
            m_components.push_back(Component(names[i], vals[i]));
}

//=============================================================================

void Composition::getComposition(
    const std::map<std::string, int>& map, double* const p_vec) const
{
    // Default values are zero
    std::fill(p_vec, p_vec+map.size(), 0.0);

    // Fill in fractions that are given
    for (int i = 0; i < size(); ++i) {
        const Component& c = m_components[i];
        std::map<std::string, int>::const_iterator iter =
            map.find(c.name);
        if (iter != map.end())
            p_vec[iter->second] = c.fraction;
        else
            throw LogicError()
                << "Component name in composition is not in the list of "
                << "possible components.";
    }
}

//=============================================================================

std::string Composition::componentsFromList(const std::string& list)
{
    // Generate list of components
    std::vector<std::string> tokens;
    String::tokenize(list, tokens, ":, \n\r\t\f");

    if (tokens.size() % 2 != 0)
        return "Composition should have the format 'name:value, name:value, ...'!";

    m_components.clear();
    for (int i = 0; i < tokens.size(); i+=2) {
        if (String::isNumeric(tokens[i+1]))
            m_components.push_back(
                Component(tokens[i], atof(tokens[i+1].c_str()))
            );
        else
            return "Composition value must be real!";
    }

    // Make sure that all the names are unique
    for (int i = 0; i < size(); ++i)
        for (int j = i+1; j < size(); ++j)
            if (m_components[i].name == m_components[j].name)
                return "Duplicate names in composition!";

    // Normalize to one
    double sum = 0.0;
    for (int i = 0; i < size(); ++i)
        sum += m_components[i].fraction;
    for (int i = 0; i < size(); ++i)
        m_components[i].fraction /= sum;

    return std::string();
}

//=============================================================================

    } // namespace Thermodynamics
} // namespace Mutation


