/**
 * @file CollisionIntegral.h
 *
 * Provides the self registering CollisionIntegral type.
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

#ifndef TRANSPORT_COLLISION_INTEGRAL_H
#define TRANSPORT_COLLISION_INTEGRAL_H

#include "Units.h"

#include <iostream>
#include <string>
#include <typeinfo>

namespace Mutation {

	// Forward declarations
	namespace Thermodynamics { class Thermodynamics; }
	namespace Utilities {
		namespace IO { class XmlElement; }
	}

	namespace Transport {

class CollisionPairNew;

/**
 * Abstract (self registering) class for all collision integral types.
 */
class CollisionIntegral
{
public:
	// Arguments for self registering constructor
	struct ARGS {
	    ARGS(
	        const Mutation::Utilities::IO::XmlElement& arg1,
	        CollisionPairNew& arg2) :
	        xml(arg1), pair(arg2) {}
	    const Mutation::Utilities::IO::XmlElement& xml;
	    CollisionPairNew& pair;
	};

	/**
	 * Constructor used in self registration.
	 */
	CollisionIntegral(ARGS args);

	/**
	 * Destructor.
	 */
	virtual ~CollisionIntegral() { };

	/// Equality operator.
	bool operator == (const CollisionIntegral& compare) const {
	    if (typeid(*this) != typeid(compare))
	        return false;
	    return isEqual(compare);
	}

	/// Equality operator.
	bool operator != (const CollisionIntegral& compare) const {
	    return !(*this == compare);
	}

	/**
	 * Returns true if this integral can be tabulated vs. temperature.  Default
	 * is to return false.
	 */
	virtual bool canTabulate() const { return false; }

	/**
	 * Returns the value of this integral at the given temperature in m^2.
	 */
	double compute(double T) {
		return m_units.convertToBase(compute_(T));
	}

	/**
	 * Gets any other parameters necessary to compute the integral which cannot
	 * be determined from the temperature alone.  Default behavior is to do
	 * nothing.
	 */
	virtual void getOtherParams(
	    const Mutation::Thermodynamics::Thermodynamics& thermo) { };

	/**
	 * Returns true if all the necessary data for this integral could be loaded.
	 */
	virtual bool loaded() const { return true; }

	/**
	 * Returns the reference for this collision integral data.
	 */
	const std::string& reference() const { return m_ref; }

	/**
	 * Returns the estimated accuracy of this collision integral in %.
	 */
	double accuracy() const { return m_acc; }

protected:

	virtual double compute_(double T) = 0;

	/**
	 * Ensures that collision integral types can be compared.
	 */
	virtual bool isEqual(const CollisionIntegral& compare) const = 0;

	/**
	 * Sets the reference for this collision integral.
	 */
	void setReference(const std::string& ref) { m_ref = ref; }

	/**
	 * Sets the estimated accuracy of this collision integral in %.
	 */
	void setAccuracy(double acc) { m_acc = acc; }

private:

	std::string m_ref;
	double      m_acc;
	Mutation::Utilities::Units m_units;
};

	} // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_COLLISION_INTEGRAL_H

