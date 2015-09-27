/**
 * @file CollisionIntegral.cpp
 *
 * Provides concrete CollisionIntegral types.
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

#include "AutoRegistration.h"
#include "CollisionIntegral.h"
#include "XMLite.h"
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;

#include <utility>
using namespace std;

namespace Mutation {
	namespace Transport {

//==============================================================================

CollisionIntegral::CollisionIntegral(CollisionIntegral::ARGS args) :
	m_ref(""),
	m_acc(0.0)
{
	XmlElement& xml = args.first;

	// Load the reference information if it exists
	xml.getAttribute("ref", m_ref, m_ref);

	// Load the estimated accuracy if that exists
	xml.getAttribute("accuracy", m_acc, m_acc);
}

//==============================================================================

/**
 * Represents a collision integral which evaluates to a constant.
 */
class ConstantColInt : public CollisionIntegral
{
public:
	ConstantColInt(CollisionIntegral::ARGS args) :
		CollisionIntegral(args)
	{
		XmlElement& xml = args.first;

		// Load the constant value
		xml.getAttribute("value", m_value,
			"A constant collision integral must provide a 'value' attribute!");
	}

	double compute(double T) { return m_value; }

private:
	double m_value;
};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<ConstantColInt, CollisionIntegral> const_ci("constant");

//==============================================================================

/**
 * Represents a collision integral which is interpolated from a given table.
 */
class TableColInt : public CollisionIntegral
{
public:
	TableColInt(CollisionIntegral::ARGS args) :
		CollisionIntegral(args)
	{
		XmlElement& xml = args.first;
		string text = xml.text();
	}

	bool canTabulate() { return true; }

	double compute(double T)
	{
		// Clip the temperature
		if (T < m_T[0])
			return m_Q[0];
		if (T > m_T.back())
			return m_Q.back();

		// Find the index
		int i = 1;
		while (m_T[i] < T && i < m_T.size()-1) i++;

		// Interpolate
		return (m_Q[i]-m_Q[i-1])*(T-m_T[i])/(m_T[i]-m_T[i-1])+m_Q[i];
	}

private:
	vector<double> m_T;
	vector<double> m_Q;
};

// Register the "table" CollisionIntegral
Config::ObjectProvider<TableColInt, CollisionIntegral> table_ci("table");

//==============================================================================

	} // namespace Transport
} // namespace Mutation
