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
#include "StringUtils.h"
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
	// Load the reference information if it exists
	args.xml.getAttribute("ref", m_ref, m_ref);

	// Load the estimated accuracy if that exists
	args.xml.getAttribute("accuracy", m_acc, m_acc);

	// Load any units that there might be
	string units; args.xml.getAttribute("units", units, string());
	if (!units.empty()) {
		// Always take last units if available
		std::vector<Units> vec = Units::split(units);
		if (vec.size() == 0) args.xml.parseError(
			"Invalid units attribute.");
		m_units = vec.back();
	}
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
		// Load the constant value
		args.xml.getAttribute("value", m_value,
			"A constant collision integral must provide a 'value' attribute!");
	}

private:

	double compute_(double T) { return m_value; }

	/**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const ConstantColInt& compare = dynamic_cast<const ConstantColInt&>(ci);
        return (m_value == compare.m_value);
    }

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
		// Get the two vectors separated by a comma
		vector<string> tokens;
		String::tokenize(args.xml.text(), tokens, ",\n\r");

		if (tokens.size() != 2) args.xml.parseError(
		    "Incorrect format for collision integral table.");

		// Parse the two rows
		fillVector(tokens[0], m_T);
		fillVector(tokens[1], m_Q);

		if (m_T.size() != m_Q.size() && m_T.size() > 1) args.xml.parseError(
		    "Table rows must be same size and greater than 1.");
	}

	// Allow tabulation of this integral type
	bool canTabulate() const { return true; }

private:

	// Interpolate from a table
	double compute_(double T)
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

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const TableColInt& compare = dynamic_cast<const TableColInt&>(ci);
        return (m_T == compare.m_T && m_Q == compare.m_Q);
    }

	/**
	 * Loads a table row from a string of numbers.
	 */
	void fillVector(const std::string& str, std::vector<double>& v) {
	    double val; istringstream ss(str);
	    while (ss >> val) v.push_back(val);
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
