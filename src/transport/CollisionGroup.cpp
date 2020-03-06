/**
 * @file CollisionGroup.cpp
 *
 * Implementation of CollisionGroup type.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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
    m_unique_vals.resize(m_size);

    // Create a list of unique integrals and a mapping from original ordering to
    // the unique list
    for (int i = 0; i < integrals.size(); ++i) {
        SharedPtr<CollisionIntegral> integral = integrals[i];
        int j;
        for (j = 0; j < m_integrals.size(); ++j)
            if (*integral == *m_integrals[j]) break;

        if (j == m_integrals.size())
            m_integrals.push_back(integral);

        m_map[i] = j;
    }

    // If not tabulating, then we are done setting up
    if (!m_tabulate)
        return;

    // Condense all the tabulatable integrals at the top of the list
    int next = 0;
    for (int i = 0; i < m_integrals.size(); ++i) {
        if (m_integrals[i]->canTabulate()) {
            // Move this integral above the non-tabulatable list
            if (next != i) {
                std::swap(m_integrals[next], m_integrals[i]);

                // Update the map
                for (int j = 0; j < m_size; ++j) {
                    if (m_map[j] == next)   m_map[j] = i;
                    else if (m_map[j] == i) m_map[j] = next;
                }
            }

            next++;
        }
    }

    // If nothing to tabulate, then we are done
    if (next == 0)
        return;

    // Generate the table
    m_table.resize(next, int((m_table_max-m_table_min)/m_table_delta)+1);

    double T = m_table_min;
    for (int j = 0; j < m_table.cols(); ++j) {
        for (int i = 0; i < m_table.rows(); ++i)
            m_table(i,j) = m_integrals[i]->compute(T);
        T += m_table_delta;
    }
}

//==============================================================================

CollisionGroup& CollisionGroup::update(
    double T, const Thermodynamics::Thermodynamics& thermo)
{
    // Compute tabulated data
    if (m_table.rows() > 0) {
        // Clip the temperature to the table bounds
        double Tc = std::max(std::min(T, m_table_max), m_table_min);

        // Compute index of temperature >= to T
        int i = std::min(
            (int)((Tc-m_table_min)/m_table_delta)+1, (int)m_table.cols()-1);
        double ratio = (Tc - m_table_min - i*m_table_delta)/m_table_delta;

        // Linearly interpolate the table
        m_unique_vals.head(m_table.rows()) =
            ratio*(m_table.col(i) - m_table.col(i-1)) + m_table.col(i);
    }

    // Compute non tabulated data
    for (int i = m_table.rows(); i < m_integrals.size(); ++i) {
        m_integrals[i]->getOtherParams(thermo);
        m_unique_vals[i] = m_integrals[i]->compute(T);
    }

    // Finally copy unique values to full vector
    for (int i = 0; i < m_size; ++i)
        m_values[i] = m_unique_vals[m_map[i]];

    return *this;
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
