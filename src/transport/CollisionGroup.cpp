/**
 * @file CollisionGroup.cpp
 *
 * Implementation of CollisionGroup type.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015 von Karman Institute for Fluid Dynamics (VKI)
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

#include "CollisionGroup.h"

#include <iostream>
using namespace std;

namespace Mutation {
    namespace Transport {

//==============================================================================

void CollisionGroup::manage(
    const std::vector< SharedPtr<CollisionIntegral> >& integrals)
{
    m_size = integrals.size();
    if (m_size == 0)
        return;

    m_map.resize(m_size);
    m_values.resize(m_size);

    // Create a list of unique integrals and a mapping from original ordering to
    // the unique list
    vector< SharedPtr<CollisionIntegral> > unique;
    for (int i = 0; i < integrals.size(); ++i) {
        SharedPtr<CollisionIntegral> integral = integrals[i];
        int j;
        for (j = 0; j < unique.size(); ++j)
            if (*integral == *unique[j]) break;

        if (j == unique.size())
            unique.push_back(integral);

        m_map[i] = j;
    }

    // If not tabulating, then we are done setting up
    if (!m_tabulate) {
        m_integrals = unique;
        return;
    }

    // Condense all the tabulatable integrals at the top of the list
    int next = 0;
    for (int i = 0; i < unique.size(); ++i) {
        if (unique[i]->canTabulate()) {
            // Move this integral above the non-tabulatable list
            if (next != i) {
                std::swap(unique[next], unique[i]);

                // Update the map
                for (int j = 0; j < m_size; ++j) {
                    if (m_map[j] == next)   m_map[j] = i;
                    else if (m_map[j] == i) m_map[j] = next;
                }
            }

            next++;
        }
    }

    cout << "managing " << m_size << " integrals (" << unique.size() << " unique)" << endl;
    cout << "  tabulating " << next << endl;
    cout << "  storing    " << unique.size()-next << endl;

    m_integrals = unique;
}

//==============================================================================

CollisionGroup& CollisionGroup::update(
    double T, const Thermodynamics::Thermodynamics& thermo)
{
    for (int i = 0; i < m_size; ++i) {
        SharedPtr<CollisionIntegral> integral = m_integrals[m_map[i]];
        integral->getOtherParams(thermo);
        m_values[i] = integral->compute(T);
    }
    return *this;
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
