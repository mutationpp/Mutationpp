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
#include "CollisionPair.h"
#include "StringUtils.h"
#include "SharedPtr.h"
#include "Thermodynamics.h"
#include "XMLite.h"
using namespace Mutation;
using namespace Mutation::Thermodynamics;
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
 * Represents a collision integral that is computed based on the expression,
 * \f[ B* = \frac{5Q^{1,2} - 4Q^{1,3}}{Q^{1,1}}, \f]
 * assumes the other data is available.
 */
class FromBstColInt : public CollisionIntegral
{
public:
    FromBstColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Determine the type of B* expression to evaluate
        string tag = args.xml.tag();
             if (tag == "Bst") m_type = BST;
        else if (tag == "Q11") m_type = Q11;
        else if (tag == "Q12") m_type = Q12;
        else if (tag == "Q13") m_type = Q13;
        else args.xml.parseError(
            "Cannot determine " + tag + " from C* expression");

        switch (m_type) {
        case BST:
            m_ci1 = args.pair.get("Q11");
            m_ci2 = args.pair.get("Q12");
            m_ci3 = args.pair.get("Q13");
            break;
        case Q11:
            m_ci1 = args.pair.get("Bst");
            m_ci2 = args.pair.get("Q12");
            m_ci3 = args.pair.get("Q13");
            break;
        case Q12:
            m_ci1 = args.pair.get("Bst");
            m_ci2 = args.pair.get("Q11");
            m_ci3 = args.pair.get("Q13");
            break;
        case Q13:
            m_ci1 = args.pair.get("Bst");
            m_ci2 = args.pair.get("Q11");
            m_ci3 = args.pair.get("Q12");
            break;
        }
    }

    /// Make sure required integrals were loaded
    bool loaded() const {
        return (m_ci1->loaded() && m_ci2->loaded() && m_ci3->loaded());
    }

    // Allow tabulation of this integral type
    bool canTabulate() const {
        return (m_ci1->canTabulate() &&
                m_ci2->canTabulate() &&
                m_ci3->canTabulate());
    }

    // Forward the getOtherParams call to the underlying integrals
    void getOtherParams(const Thermodynamics::Thermodynamics& thermo) {
        m_ci1->getOtherParams(thermo);
        m_ci2->getOtherParams(thermo);
        m_ci3->getOtherParams(thermo);
    }

private:

    enum BstType {
        BST, Q11, Q12, Q13
    };

    // Evaluate the B* expression
    double compute_(double T) {
        switch(m_type) {
        case BST: return // (5Q12 - 4Q13)/Q11
            (5.*m_ci2->compute(T)-4.*m_ci3->compute(T))/m_ci1->compute(T);
        case Q11: return // (5Q12 - 4Q13)/B*
            (5.*m_ci2->compute(T)-4.*m_ci3->compute(T))/m_ci1->compute(T);
        case Q12: return // (B*Q11 + 4Q13)/5
            (m_ci1->compute(T)*m_ci2->compute(T)+4.*m_ci3->compute(T))/5.;
        case Q13: return // (5*Q12 - B*Q11)/4
            (5.*m_ci3->compute(T)-m_ci1->compute(T)*m_ci2->compute(T))/4.;
        }
    }

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const FromBstColInt& compare = dynamic_cast<const FromBstColInt&>(ci);
        return (*m_ci1 == *(compare.m_ci1) &&
                *m_ci2 == *(compare.m_ci2) &&
                *m_ci3 == *(compare.m_ci3));
    }

private:

    BstType m_type;

    SharedPtr<CollisionIntegral> m_ci1;
    SharedPtr<CollisionIntegral> m_ci2;
    SharedPtr<CollisionIntegral> m_ci3;
};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<FromBstColInt, CollisionIntegral> bst_ci("from B*");

//==============================================================================

/**
 * Represents a collision integral that is computed based on the expression,
 * \f[ C* = \frac{Q^{1,2}}{Q^{1,1}}, \f]
 * assumes the other data is available.
 */
class FromCstColInt : public CollisionIntegral
{
public:
    FromCstColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Determine the type of B* expression to evaluate
        string tag = args.xml.tag();
             if (tag == "Cst") m_type = CST;
        else if (tag == "Q11") m_type = Q11;
        else if (tag == "Q12") m_type = Q12;
        else args.xml.parseError(
            "Cannot determine " + tag + " from C* expression");

        switch (m_type) {
        case CST:
            m_ci1 = args.pair.get("Q11");
            m_ci2 = args.pair.get("Q12");
            break;
        case Q11:
            m_ci1 = args.pair.get("Cst");
            m_ci2 = args.pair.get("Q12");
            break;
        case Q12:
            m_ci1 = args.pair.get("Cst");
            m_ci2 = args.pair.get("Q11");
            break;
        }
    }

    /// Make sure required integrals were loaded
    bool loaded() const { return (m_ci1->loaded() && m_ci2->loaded()); }

    // Allow tabulation of this integral type
    bool canTabulate() const {
        return (m_ci1->canTabulate() && m_ci2->canTabulate());
    }

    // Forward the getOtherParams call to the underlying integrals
    void getOtherParams(const Thermodynamics::Thermodynamics& thermo) {
        m_ci1->getOtherParams(thermo);
        m_ci2->getOtherParams(thermo);
    }

private:

    enum CstType {
        CST, Q11, Q12
    };

    // Evaluate the B* expression
    double compute_(double T) {
        switch(m_type) {
        case CST: return // Q12/Q11
            (m_ci2->compute(T)/m_ci1->compute(T));
        case Q11: return // Q12/C*
            (m_ci2->compute(T)/m_ci1->compute(T));
        case Q12: return // C* Q11
            (m_ci1->compute(T)*m_ci2->compute(T));
        }
    }

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const FromCstColInt& compare = dynamic_cast<const FromCstColInt&>(ci);
        return (*m_ci1 == *(compare.m_ci1) && *m_ci2 == *(compare.m_ci2));
    }

private:

    CstType m_type;

    SharedPtr<CollisionIntegral> m_ci1;
    SharedPtr<CollisionIntegral> m_ci2;
};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<FromCstColInt, CollisionIntegral> cst_ci("from C*");

//==============================================================================

/**
 * Represents a collision integral which is taken as the ratio of another
 * integral.
 */
class RatioColInt : public CollisionIntegral
{
public:
    RatioColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Load the ratio
        args.xml.getAttribute("ratio", m_ratio,
            "A ratio collision integral must provide a 'ratio' attribute.");

        // Load the integral
        string integral;
        args.xml.getAttribute("integral", integral,
            "A ratio collision integral must provide a 'integral' attribute.");
        m_integral = args.pair.get(integral);
    }

    // Allow tabulation of this integral type
    bool canTabulate() const { return m_integral->canTabulate(); }

    // Forward the getOtherParams call to the underlying integral
    void getOtherParams(const Thermodynamics::Thermodynamics& thermo) {
        m_integral->getOtherParams(thermo);
    }

private:

    double compute_(double T) { return m_ratio * m_integral->compute(T); }

    /**
     * Returns true if the ratio and integral are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const RatioColInt& compare = dynamic_cast<const RatioColInt&>(ci);
        return (m_ratio == compare.m_ratio &&
            *m_integral == *(compare.m_integral));
    }

private:

    double m_ratio;
    SharedPtr<CollisionIntegral> m_integral;
};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<RatioColInt, CollisionIntegral> ratio_ci("ratio");

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
