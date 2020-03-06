/**
 * @file CollisionIntegral.cpp
 *
 * Provides concrete CollisionIntegral types.
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

#include "AutoRegistration.h"
#include "CollisionIntegral.h"
#include "CollisionPair.h"
#include "Constants.h"
#include "Interpolators.h"
#include "StringUtils.h"
#include "SharedPtr.h"
#include "Thermodynamics.h"
#include "XMLite.h"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Mutation::Numerics;

#include <utility>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cassert>
#include <sstream>

using namespace std;

namespace Mutation {
	namespace Transport {

//==============================================================================

SharedPtr<CollisionIntegral> CollisionIntegral::load(
    CollisionIntegral::ARGS args)
{
    string type;
    args.xml.getAttribute("type", type, "Integral type must be specified.");

    CollisionIntegral* p_ci;
    try {
        p_ci = Config::Factory<CollisionIntegral>::create(type, args);
    } catch (Error& e) {
        e << "\nWas trying to load a collision integral.";
        throw;
    }

    args.xml.parseCheck(p_ci != NULL,
        "Invalid collision integral type '" + type + "' for " + args.kind +
        "_" + args.pair.name() + ".");

    args.xml.parseCheck(p_ci->loaded(),
        "Could not find data necessary to load " + args.kind +
        "_" + args.pair.name() + ".");

    return SharedPtr<CollisionIntegral>(p_ci);
}

//==============================================================================

CollisionIntegral::CollisionIntegral(CollisionIntegral::ARGS args) :
	m_ref(""),
	m_acc(0.0),
	m_fac(1.0)
{
	// Load the reference information if it exists
	args.xml.getAttribute("ref", m_ref, m_ref);

	// Load the estimated accuracy if that exists
	args.xml.getAttribute("accuracy", m_acc, m_acc);

	// Load any units that there might be
	string units; args.xml.getAttribute("units", units, units);
	if (!units.empty()) {
		// Always take last units if available
		std::vector<Units> vec = Units::split(units);
		if (vec.size() == 0) args.xml.parseError(
			"Invalid units attribute.");
		m_units = vec.back();
	}

	// Check if we should multiply by pi
	bool multpi = false; args.xml.getAttribute("multpi", multpi, multpi);
	m_fac = (multpi ? PI : 1.0);
}

//==============================================================================

double CollisionIntegral::loadSpeciesParameter(
    Mutation::Utilities::IO::XmlElement& root, const std::string& parameter,
    const std::string& species, std::string units, double def)
{
    // First look for the parameter table
    XmlElement::const_iterator db_iter =
        root.findTag(parameter);
    if (db_iter == root.end())
        return def;

    // Load default units if available
    db_iter->getAttribute("units", units, units);

    // Check if the species name has an alias
    std::string name = speciesAlias(root, species);

    // Look for species name
    XmlElement::const_iterator sp_iter =
        db_iter->findTagWithAttribute("species", "name", name);
    if (sp_iter == db_iter->end())
        return def;

    // Found it, get value
    double value;
    sp_iter->getAttribute("value", value, "must have a 'value' attribute.");

    // Check if units are different then the default
    sp_iter->getAttribute("units", units, units);

    // Return value in base units
    if (units.empty())
        return value;

    return Units(units).convertToBase(value);
}

//==============================================================================

std::string CollisionIntegral::speciesAlias(
	Mutation::Utilities::IO::XmlElement& root, const std::string& species)
{
    // First look for the aliases list
    XmlElement::const_iterator db_iter =
        root.findTag("aliases");
    if (db_iter == root.end())
        return species;
    
    // Look for species name
    XmlElement::const_iterator sp_iter =
        db_iter->findTagWithAttribute("species", "name", species);
    if (sp_iter == db_iter->end())
        return species;
    
    // Get the alias name
    std::string alias;
    sp_iter->getAttribute("alias", alias, "must specify a species alias");
    return alias;
}

//==============================================================================

/**
 * Represents a constant value collision integral which writes a warning that
 * the requested integral is not found in the database.
 */
class WarningColInt : public CollisionIntegral
{
public:
    WarningColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Load the constant value
        args.xml.getAttribute("value", m_value,
            "A warning collision integral must provide a 'value' attribute!");

        // Write a warning
        cout << "Warning: missing collision integral " << args.xml.tag() << "_("
             << args.pair.sp1Name() << "," << args.pair.sp2Name()
             << ").  Using a constant value of " << m_value << "." << endl;
    }

private:

    double compute_(double T) { return m_value; }

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const WarningColInt& compare = dynamic_cast<const WarningColInt&>(ci);
        return (m_value == compare.m_value);
    }

private:

    double m_value;
};

// Register the "warning" CollisionIntegral
Config::ObjectProvider<WarningColInt, CollisionIntegral> warn_ci("warning");

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
 * @brief A collision integral evaluated as an exponential polynomial.
 * 
 * The general exponential polynomial is given as 
 * \f[
 *     f(T) = \exp( \sum_{i=0}^{n-1} A_i (\ln T)^i ),
 * \f]
 * where \f$A_i\f$ are the polynomial coefficients.
 */
class ExpPolyColInt : public CollisionIntegral
{
public:
    ExpPolyColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Load polynomial coefficients
        std::stringstream ss(args.xml.text());

        std::copy(
            std::istream_iterator<double>(ss),
            std::istream_iterator<double>(),
            std::back_inserter(m_params));
    }

private:

    double compute_(double T) {
        double lnT = std::log(T);
        double val = m_params[0];
        for (int i = 1; i < m_params.size(); ++i)
            val = val*lnT + m_params[i];
        return std::exp(val);
    }

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const ExpPolyColInt& compare = dynamic_cast<const ExpPolyColInt&>(ci);
        return (m_params == compare.m_params);
    }

private:

    std::vector<double> m_params;

};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<ExpPolyColInt, CollisionIntegral> exppoly_ci("exp-poly");

//==============================================================================


/**
 * Represents a collision integral that is computed based on the expression,
 * \f[ A* = \frac{Q^{2,2}}{Q^{1,1}}, \f]
 * assumes the other data is available.
 */
class FromAstColInt : public CollisionIntegral
{
public:
    FromAstColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Determine the type of B* expression to evaluate
        string tag = args.xml.tag();
             if (tag == "Ast") m_type = AST;
        else if (tag == "Q11") m_type = Q11;
        else if (tag == "Q22") m_type = Q22;
        else args.xml.parseError(
            "Cannot determine " + tag + " from A* expression");

        switch (m_type) {
        case AST:
            m_ci1 = args.pair.get("Q11");
            m_ci2 = args.pair.get("Q22");
            break;
        case Q11:
            m_ci1 = args.pair.get("Ast");
            m_ci2 = args.pair.get("Q22");
            break;
        case Q22:
            m_ci1 = args.pair.get("Ast");
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

    enum AstType {
        AST, Q11, Q22
    };

    // Evaluate the A* expression
    double compute_(double T) {
        switch(m_type) {
        case AST: return // Q22/Q11
            (m_ci2->compute(T)/m_ci1->compute(T));
        case Q11: return // Q22/A*
            (m_ci2->compute(T)/m_ci1->compute(T));
        case Q22: return // A* Q11
            (m_ci1->compute(T)*m_ci2->compute(T));
        }
    }

    /**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const FromAstColInt& compare = dynamic_cast<const FromAstColInt&>(ci);
        return (*m_ci1 == *(compare.m_ci1) && *m_ci2 == *(compare.m_ci2));
    }

private:

    AstType m_type;

    SharedPtr<CollisionIntegral> m_ci1;
    SharedPtr<CollisionIntegral> m_ci2;
};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<FromAstColInt, CollisionIntegral> ast_ci("from A*");

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
        // Determine the type of C* expression to evaluate
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

    // Evaluate the A* expression
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
 * Represents a collision integral computed using Murphy's phenomenological
 * approximation to combine elastic and charge-exchange integrals.  The
 * approximation is given as
 * \f[
 *    f(T) = \sqrt(Q_1(T)^2 + Q_2(T)^2)
 * \f]
 * where \f$\sqrt(Q_1)\f$ and \f$\sqrt(Q_2)\f$ are the elastic and
 * charge-exchange integrals, respectively.
 */
class MurphyColInt : public CollisionIntegral
{
public:
    MurphyColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        // Load the two integrals
        m_Q1 = getIntegral(args, "Q1");
        m_Q2 = getIntegral(args, "Q2");
    }

    // Allow tabulation of this integral if the two underlying integrals can be
    // tabulated
    bool canTabulate() const {
        return m_Q1->canTabulate() && m_Q2->canTabulate();
    }

    /**
     * Returns true if the ratio and integral are the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const MurphyColInt& compare = dynamic_cast<const MurphyColInt&>(ci);
        return (*m_Q1 == *(compare.m_Q1) &&
                *m_Q2 == *(compare.m_Q2));
    }


private:

    double compute_(double T)
    {
        double Q1 = m_Q1->compute(T);
        double Q2 = m_Q2->compute(T);
        return std::sqrt(Q1*Q1 + Q2*Q2);
    }

    SharedPtr<CollisionIntegral> getIntegral(
            CollisionIntegral::ARGS args, const std::string& tag)
    {
        XmlElement::const_iterator Q = args.xml.findTag(tag);
        if (Q == args.xml.end())
            args.xml.parseError(
                "Murphy collision integral must have " + tag + " child node.");

        return CollisionIntegral::load(
            CollisionIntegral::ARGS(*Q, args.pair, args.kind));
    }

private:

    SharedPtr<CollisionIntegral> m_Q1;
    SharedPtr<CollisionIntegral> m_Q2;
};

// Register the "constant" CollisionIntegral
Config::ObjectProvider<MurphyColInt, CollisionIntegral> murphy_ci("Murphy");

//==============================================================================

/**
 * @brief Represents a collision integral which is interpolated from a table.
 * 
 * The default interpolation scheme is linear, however this can be changed
 * by setting the interpolator attribute of the integral.  In addition, the
 * default behavior is to clip the integral value at the bounds of the table
 * (min and max temperature) to the values at those temperatures.  This too can
 * be turned off by setting the the clip attribute to false or no.
 */
class TableColInt : public CollisionIntegral
{
public:
	TableColInt(CollisionIntegral::ARGS args) :
		CollisionIntegral(args),
		m_interpolator_type("Linear"), m_clip(true)
	{
		// First load global options
	    loadGlobalOptions(args.xml.document()->root());

	    // Check for the same options here
	    // clip
        args.xml.getAttribute("clip", m_clip, m_clip);

        // interpolator
        args.xml.getAttribute(
            "interpolator", m_interpolator_type, m_interpolator_type);

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

		// Load the interpolator
		try {
            mp_interpolator = SharedPtr< Interpolator<double> > (
                Config::Factory<Interpolator<double> >::create(
                    m_interpolator_type,
                    Interpolator<double>::ARGS(&m_T[0], &m_Q[0], m_T.size())
                ));
		} catch (Error& e) {
		    e << "\nWas trying to load a collision integral interpolation "
		      << "scheme.";
		    throw;
		}
	}

	// Allow tabulation of this integral type
	bool canTabulate() const { return true; }

private:

	// Interpolate from a table
	double compute_(double T)
	{
		// Clip the temperature if requested
	    if (m_clip) {
            if (T < m_T[0])
                return m_Q[0];
            if (T > m_T.back())
                return m_Q.back();
	    }

		return (*mp_interpolator)(T);
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

	void loadGlobalOptions(XmlElement& root)
	{
	    // Look for global-options element
        XmlElement::const_iterator options = root.findTag("global-options");
        if (options == root.end()) return;

        // Look for table element
        XmlElement::const_iterator table =
            options->findTagWithAttribute("integral", "type", "table");
        if (table == options->end()) return;

        // clip
        table->getAttribute("clip", m_clip, m_clip);

        // interpolator
        table->getAttribute(
            "interpolator", m_interpolator_type, m_interpolator_type);
	}

private:

	vector<double> m_T;
	vector<double> m_Q;

	SharedPtr< Interpolator<double> > mp_interpolator;
    string m_interpolator_type;
	bool m_clip;

};

// Register the "table" CollisionIntegral
Config::ObjectProvider<TableColInt, CollisionIntegral> table_ci("table");

//==============================================================================

	} // namespace Transport
} // namespace Mutation
