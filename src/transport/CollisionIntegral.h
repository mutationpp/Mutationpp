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

#include <string>

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
	typedef const std::pair<Utilities::IO::XmlElement&, CollisionPair&>& ARGS;

	/**
	 * Constructor used in self registration.
	 */
	CollisionIntegral(ARGS args);

	/**
	 * Destructor.
	 */
	virtual ~CollisionIntegral() { };

	/**
	 * Returns true if this integral can be tabulated vs. temperature.
	 */
	virtual bool canTabulate() { return false; }

	/**
	 * Returns the value of this integral at the given temperature in m^2.
	 */
	virtual double compute(double T) = 0;

	/**
	 * Gets any other parameters necessary to compute the integral which cannot
	 * be determined from the temperature alone.
	 */
	virtual void getOtherParams(const Thermodynamics::Thermodynamics& thermo) { };

	/**
	 * Returns true if all the necessary data for this integral could be loaded.
	 */
	virtual bool loaded() { return true; }

	/**
	 * Returns the reference for this collision integral data.
	 */
	const std::string& reference() const { return m_ref; }

	/**
	 * Returns the estimated accuracy of this collision integral in %.
	 */
	double accuracy() const { return m_acc; }

protected:

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
};

	} // namespace Transport
} // namespace Mutation


