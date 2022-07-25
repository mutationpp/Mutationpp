/**
 * @file CollisionIntegral.h
 *
 * Provides the self registering CollisionIntegral type.
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

#ifndef TRANSPORT_COLLISION_INTEGRAL_H
#define TRANSPORT_COLLISION_INTEGRAL_H

#include "Units.h"
#include "SharedPtr.h"

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

class CollisionPair;

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
	        CollisionPair& arg2, const std::string& arg3) :
	        xml(arg1), pair(arg2), kind(arg3) {}
	    const Mutation::Utilities::IO::XmlElement& xml;
	    CollisionPair& pair;
	    std::string kind;
	};

	/// Returns name of this type.
	static std::string typeName() { return "CollisionIntegral"; }

	/**
	 * Constructor used in self registration.
	 */
	CollisionIntegral(ARGS args);

	/**
	 * Destructor.
	 */
	virtual ~CollisionIntegral() { };

	/// Equality operator.
	virtual bool operator == (const CollisionIntegral& compare) const {
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
		return m_fac*m_units.convertToBase(compute_(T));
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

	/**
	 * Tries to load the integral.
	 */
	static SharedPtr<CollisionIntegral> load(ARGS args);

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

	/**
	 * Sets the units of the integral.
	 */
	void setUnits(const Mutation::Utilities::Units& units) { m_units = units; }

	/**
	 * Sets the multiplication factor of this integral.
	 */
	void setFactor(double fac) { m_fac = fac; }

	/**
	 * Loads pure species parameter from the collision database with the given
	 * name for the given species in base units.  If the parameter is not
	 * found, then the default value is returned.
	 */
	double loadSpeciesParameter(
	    Mutation::Utilities::IO::XmlElement& root,
	    const std::string& parameter,
	    const std::string& species,
	    const std::string units = "",
	    double def = -1.0);
	
	/**
	 * Returns the alias of this species in the collision integral database.
	 * If there is no alias found, then the name is returned.
	 * 
	 * @param root Root XMLElement representing the collision integral
	 * database.
	 * @param species The species.
	 * @return std::string The alias name.
	 */
	std::string speciesAlias(
		Mutation::Utilities::IO::XmlElement& root,
		const std::string& species);

private:

	std::string m_ref;
	double m_acc;
	double m_fac;
	Mutation::Utilities::Units m_units;
};

	} // namespace Transport
} // namespace Mutation

#endif // TRANSPORT_COLLISION_INTEGRAL_H

