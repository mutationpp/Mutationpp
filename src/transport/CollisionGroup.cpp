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

namespace Mutation {
    namespace Transport {

//==============================================================================

void CollisionGroup::manage(
    const std::vector< SharedPtr<CollisionIntegral> >& integrals)
{
    m_integrals = integrals;
    m_size = m_integrals.size();
    m_values.resize(m_size);

}

//==============================================================================

CollisionGroup& CollisionGroup::update(
    double T, const Thermodynamics::Thermodynamics& thermo)
{
    for (int i = 0; i < m_integrals.size(); ++i) {
        m_integrals[i]->getOtherParams(thermo);
        m_values[i] = m_integrals[i]->compute(T);
    }
    return *this;
}

//==============================================================================

    } // namespace Transport
} // namespace Mutation
