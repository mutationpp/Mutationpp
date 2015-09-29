/**
 * @file CoulombIntegrals.cpp
 *
 * Provides Coulomb collision integral types.
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
#include "Constants.h"
#include "CollisionIntegral.h"
#include "CollisionPair.h"
#include "Thermodynamics.h"
#include "XMLite.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities::Config;

#include <Eigen/Dense>
using namespace Eigen;

#include <iostream>
#include <string>
using namespace std;

namespace Mutation {
	namespace Transport {

enum CoulombType {
	Q11_ATT = 0,
	Q11_REP,
	Q22_ATT,
	Q22_REP,
	BST_ATT,
	BST_REP,
	CST_ATT,
	CST_REP,
	N_COULOMB_TYPES // this must stay last
};

/**
 * Evaluates the Debye-Huckle collision integrals all at once because this is
 * faster then doing them individually all the time.
 */
class DebyeHuckleEvaluator
{
public:
	DebyeHuckleEvaluator() :
		m_lambda(0.0), m_last_T(-1.0), m_last_lambda(-1.0)
	{
		// Coefficients used to compute the exponential polynomials
		m_coeffs <<
			-7.927047e-1,+5.886772e-1,-8.760713e-2,+9.301813e-3,-4.131521e-4, // Q11_ATT
			-1.398075e+0,+8.048207e-1,-9.480165e-2,+5.281218e-3,-8.265206e-5, // Q11_REP
			-8.114574e-1,+7.041926e-1,-1.221972e-1,+1.323471e-2,-5.699409e-4, // Q22_ATT
			-1.108917e+0,+7.746086e-1,-1.016815e-1,+6.760088e-3,-1.562221e-4, // Q22_REP
			+2.861634e-1,-4.885811e-2,+1.337277e-3,+5.668631e-4,-4.496418e-5, // BST_ATT
			+3.260697e-1,-2.374246e-2,-9.698775e-3,+1.779920e-3,-8.393875e-5, // BST_REP
			-6.467276e-1,-1.073724e-1, 1.796216e-2,-1.865304e-3,+8.013450e-5, // CST_ATT
			-5.007713e-1,-1.063984e-1,-1.147815e-3,+1.835593e-3,-1.183034e-4; // CST_REP
	}

	double operator () (double T, CoulombType type)
	{
		// Check if we need to update the values
		if (std::abs(T-m_last_T) + std::abs(m_lambda-m_last_lambda) > 1.0e-10) {
			// Average closest impact parameter
			const double b = QE * QE / (8.0 * PI * EPS0 * KB) / T;

			// Reduced temperature vector
			m_tstvec[0] = 1.0;
			m_tstvec[1] = std::log(std::max(0.5 * m_lambda / b, 0.1));
			m_tstvec[2] = m_tstvec[1]*m_tstvec[1];
			m_tstvec[3] = m_tstvec[2]*m_tstvec[2];
			m_tstvec[4] = m_tstvec[3]*m_tstvec[3];

			// Evaluate the polynomials with mat/vec multiply
			m_values.matrix() = m_coeffs * m_tstvec;

			// Compute exponentials
			m_values = m_values.exp();

			// Dimensionalize from reduced form
			m_values.head(BST_ATT) *=
				PI*m_lambda*m_lambda/(m_tstvec[2]*m_tstvec[2]);

			// Save last updated state
			m_last_T = T;
			m_last_lambda = m_lambda;
		}

		return m_values[type];
	}

	/**
	 * Sets the Debye length given electron temperature and number density.
	 */
	void setDebyeLength(double Te, double ne) {
		assert(Te >= 0.0);
		assert(ne >= 0.0);
		m_lambda = (ne > 0.0 ? std::sqrt(EPS0*KB*Te/(2.0*ne*QE*QE)) : 0.0);
	}

private:

	double m_lambda;
	double m_last_T;
	double m_last_lambda;

	Matrix<double, N_COULOMB_TYPES, 5> m_coeffs;
	Array< double, N_COULOMB_TYPES, 1> m_values;
	Matrix<double, 5, 1>               m_tstvec;
};

/**
 * Represents a collision integral computed with the Debye-Huckle potential
 * (Coulomb potential screened by the Debye length).
 */
class DebyeHuckleColInt : public CollisionIntegral
{
public:

	/**
	 * Self registering onstructor. Determines which integral to evaluate based
	 * on the collision type.
	 */
	DebyeHuckleColInt(CollisionIntegral::ARGS args) :
		CollisionIntegral(args)
	{
		assert(args.pair.type() == ATTRACTIVE || args.pair.type() == REPULSIVE);

		// Set the type of integral we need
		string tag = args.xml.tag();
		CollisionType type = args.pair.type();

		if (tag == "Q11")
			m_type = (type == ATTRACTIVE ? Q11_ATT : Q11_REP);
		else if (tag == "Q22")
			m_type = (type == ATTRACTIVE ? Q22_ATT : Q22_REP);
		else if (tag == "Bst")
			m_type = (type == ATTRACTIVE ? BST_ATT : BST_REP);
		else if (tag == "Cst")
			m_type = (type == ATTRACTIVE ? CST_ATT : CST_REP);
		else
			args.xml.parseError(
				"Invalid collision integral for Debye-Huckle integral.");
	}

	// Set the Debye length
	void getOtherParams(const class Thermodynamics& thermo) {
		sm_evaluator.setDebyeLength(
			thermo.Te(),
			thermo.hasElectrons() ? thermo.numberDensity()*thermo.X()[0] : 0.0);
	};

private:

	double compute_(double T) {	return sm_evaluator(T, m_type);	}

	/**
     * Returns true if the constant value is the same.
     */
    bool isEqual(const CollisionIntegral& ci) const {
        const DebyeHuckleColInt& compare =
        	dynamic_cast<const DebyeHuckleColInt&>(ci);
        return (m_type == compare.m_type);
    }

private:

    CoulombType m_type;
    static DebyeHuckleEvaluator sm_evaluator;

}; // class CoulombColInt

// Initialization of the DebyeHuckleEvaluator
DebyeHuckleEvaluator DebyeHuckleColInt::sm_evaluator;

// Register the "Debye-Huckle" CollisionIntegral
ObjectProvider<DebyeHuckleColInt, CollisionIntegral> DebyeHuckle_ci("Debye-Huckle");

	} // namespace Transport
} // namespace Mutation
