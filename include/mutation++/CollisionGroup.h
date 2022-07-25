/**
 * @file CollisionGroup.h
 *
 * Provides CollisionGroup type.
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

#ifndef TRANSPORT_COLLISION_GROUP_H
#define TRANSPORT_COLLISION_GROUP_H

#include "CollisionIntegral.h"
#include "SharedPtr.h"

#include <Eigen/Dense>

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
     * Constructs an empty CollisionGroup with specified tabulation parameters.
     * @param tabulate - whether or not to tabulate collision integrals
     * @param min      - minimum temperature for tabulation
     * @param max      - maximum temperature for tabulation
     * @param delta    - temperature spacing in the table
     */
    CollisionGroup(
        bool tabulate = true,
        double min = 300.0, double max = 20000.0, double delta = 100.0) :
        m_tabulate(tabulate),
        m_size(0),
        m_table_min(min), m_table_max(max), m_table_delta(delta)
    { }

    /**
     * Sets the collision integrals that are managed by this group.
     */
    template <typename Iterator, typename Provider, typename Params>
    void manage(
        Iterator first, Iterator last,
        SharedPtr<CollisionIntegral> (Provider::*f)(const Params&),
        const Params& params)
    {
        std::vector< SharedPtr<CollisionIntegral> > integrals;
        while (first != last)
            integrals.push_back(((*first++).*f)(params));
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

    /**
     * Access the value of the ith collision integral.
     */
    const double& operator () (int i) const { return m_values[i]; }

    /**
     * Access the entire array of collision integral values.
     */
    const Eigen::ArrayXd& array() const { return m_values; }

private:

    /// Whether or not to tabulate integrals managed by this group
    bool m_tabulate;

    /// Number of integrals managed by this group
    int  m_size;

    /// vector of non-tabulated integrals
    std::vector< SharedPtr<CollisionIntegral> > m_integrals;

    /// Internal vector of computed collision integral values
    Eigen::ArrayXd   m_values;
    Eigen::ArrayXd   m_unique_vals;
    std::vector<int> m_map;

    /// Table of tabulated integrals versus temperature
    double m_table_min;
    double m_table_max;
    double m_table_delta;
    Eigen::ArrayXXd m_table;
};

	} // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_COLLISION_GROUP_H
