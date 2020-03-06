/**
 * @file CoulombIntegral.h
 *
 * Exposes the DebyeHuckleEvaluator type in case a user is interested in using
 * this class as a stand-alone object.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2016-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef TRANSPORT_COULOMB_INTEGRAL_H
#define TRANSPORT_COULOMB_INTEGRAL_H

#include <Eigen/Dense>

namespace Mutation {
    namespace Transport {

/**
 * List of collision integral types that can be computed using predefined
 * collision integrals computed from the screened Coulomb potential.
 */
enum CoulombType {
    Q11_ATT, Q11_REP,
    Q22_ATT, Q22_REP,
    Q14_ATT, Q14_REP,
    Q15_ATT, Q15_REP,
    Q24_ATT, Q24_REP,
    BST_ATT, BST_REP,
    CST_ATT, CST_REP,
    EST_ATT, EST_REP,
    N_TABLE_TYPES, // All collision integrals in the table must come before this
    Q12_ATT, Q12_REP,
    Q13_ATT, Q13_REP,
    Q23_ATT, Q23_REP,
    AST_ATT, AST_REP
};

// Number of temperature points in the collision integral table
#define N_TST_POINTS 26

/**
 * Evaluates the Debye-Huckle collision integrals all at once because this is
 * faster then doing them individually all the time.
 */
class DebyeHuckleEvaluator
{
public:
    /**
     * Constructor.
     */
    DebyeHuckleEvaluator() :
        m_lambda(0.0), m_last_T(-1.0), m_last_lambda(-1.0)
    { }

    /**
     * Computes the Coulomb integral type at the temperature T.
     */
    double operator () (double T, CoulombType type);

    /**
     * Sets the Debye length given electron temperature and number density.
     */
    void setDebyeLength(double Te, double ne);

private:

    /**
     * Interpolates the table based on reduced temperature.
     */
    void interpolate(double Tst);

private:

    double m_lambda;
    double m_last_T;
    double m_last_lambda;

    static Eigen::Array<double, N_TST_POINTS, N_TABLE_TYPES, Eigen::ColMajor> sm_table;
    static Eigen::Array<double, N_TST_POINTS, 1> sm_tstvec;
    Eigen::Array<double, N_TABLE_TYPES, 1> m_values;
};


    } // namespace Transport
} // namespace Mutation


#endif // TRANSPORT_COULOMB_INTEGRAL_H
