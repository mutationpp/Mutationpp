/**
 * @file Composition.h
 *
 * @brief Declaration of the Thermodynamics::Composition class.
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

#ifndef THERMO_COMPOSITION_H
#define THERMO_COMPOSITION_H

#include <string>
#include <vector>
#include <map>
#include <cassert>

#include "XMLite.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * Represents a single, named, composition with (component name,fraction) pairs.
 * Provides a convenient object for representing elemental or species
 * compositions with consistent error checking and normalization.
 */
class Composition
{
public:
    enum Type {
        MASS,
        MOLE
    };

    struct Component {
        Component(const std::string& n, double f)
            : name(n), fraction(f)
        { }
        std::string name;
        double fraction;
    };

    // Constructors
    Composition(const char* const list);
    Composition(
        const std::string& name, const std::string& list, Type type = MOLE);
    Composition(const Mutation::Utilities::IO::XmlElement& element);
    Composition(
        const std::vector<std::string>& names,
        const double* const vals, Type type);

    // Property getters
    const std::string& name() const { return m_name; }
    Type type() const { return m_type; }
    const Component& operator [] (int i) const {
        assert(i >= 0);
        assert(i < size());
        return m_components[i];
    }
    void getComposition(
        const std::map<std::string, int>& map, double* const p_vec) const;
    int size() const { return m_components.size(); }


private:

    std::string componentsFromList(const std::string& list);

private:

    std::string m_name;
    Type m_type;
    std::vector<Component> m_components;
};

    } // namespace Thermodynamics
} // namespace Mutation

#endif // THERMO_COMPOSITION_H
