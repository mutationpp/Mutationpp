/**
 * @file CollisionGroup.h
 *
 * Provides CollisionGroup type.
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

#ifndef TRANSPORT_COLLISION_GROUP_H
#define TRANSPORT_COLLISION_GROUP_H

#include "CollisionIntegral.h"
#include "SharedPtr.h"

#include <vector>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}

namespace Mutation {
	namespace Transport {

/**
 * Manages a collection of collision integrals which are evaluated at the same
 * temperature and conditions.
 */
class CollisionGroup
{
public:

    /**
     * Empty constructor.
     */
    CollisionGroup() : m_tabulate(true), m_size(0) { }

    /**
     * Sets the collision integrals that are managed by this group.
     */
    template <typename Iterator, typename Provider>
    void manage(
        Iterator first, Iterator last,
        SharedPtr<CollisionIntegral> (Provider::*f)() const)
    {
        std::vector< SharedPtr<CollisionIntegral> > integrals;
        while (first != last)
            integrals.push_back(((*first++).*f)());
        manage(integrals);
    }

    /**
     * Sets the collision integrals that are managed by this group.
     */
    void manage(const std::vector< SharedPtr<CollisionIntegral> >& integrals);

    /**
     * Updates the collision integral values for this collision group using the
     * given temperature.  Returns a reference to itself.
     */
    CollisionGroup& update(
        double T, const Thermodynamics::Thermodynamics& thermo);

    /**
     * Number of integrals managed by this group.
     */
    int size() const { return m_size; }

    /**
     * Access the value of the ith collision integral.
     */
    const double& operator [] (int i) const { return m_values[i]; }

private:

    bool m_tabulate;
    int  m_size;

    /// vector of non-tabulated integrals
    std::vector< SharedPtr<CollisionIntegral> > m_integrals;
    std::vector<double> m_values;
};

	} // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_COLLISION_GROUP_H
