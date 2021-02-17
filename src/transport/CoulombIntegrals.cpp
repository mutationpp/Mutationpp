/**
 * @file CoulombIntegrals.cpp
 *
 * Provides Coulomb collision integral types.
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
#include "Constants.h"
#include "CollisionIntegral.h"
#include "CollisionPair.h"
#include "CoulombIntegrals.h"
#include "Thermodynamics.h"
#include "XMLite.h"
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities::Config;

#include <Eigen/Dense>
using namespace Eigen;

#include <cassert>
#include <string>
using namespace std;

namespace Mutation {
    namespace Transport {


double DebyeHuckleEvaluator::operator() (double T, CoulombType type)
{
    // Check input data
    assert(T >= 0.0);
    assert(type != N_TABLE_TYPES);

    // Average closest impact parameter
    const double b = QE*QE / (8.0*PI*EPS0*KB*T);

    // Clip maximum lambda
    m_lambda = std::min(m_lambda, 2.0*sm_tstvec[N_TST_POINTS-1]*b);

    // Check if we need to update the values
    if (std::abs(T-m_last_T) + std::abs(m_lambda-m_last_lambda) > 1.0e-10) {
        // Reduced temperature
        const double Tst = std::max(0.5*m_lambda/b, sm_tstvec[0]);

        // Interpolate the table using reduced temperature
        interpolate(Tst);

        // Convert from reduced form
        m_values.head(BST_ATT) *= PI*m_lambda*m_lambda/(Tst*Tst);

        // Save last updated state
        m_last_T = T;
        m_last_lambda = m_lambda;
    }

    switch (type) {
    // Derived data from the tabulated ones
    case Q12_ATT: return m_values[CST_ATT]*m_values[Q11_ATT];
    case Q12_REP: return m_values[CST_REP]*m_values[Q11_REP];
    case Q13_ATT: return m_values[Q11_ATT]*(1.25*m_values[CST_ATT]-0.25*m_values[BST_ATT]);
    case Q13_REP: return m_values[Q11_REP]*(1.25*m_values[CST_REP]-0.25*m_values[BST_REP]);
    case Q23_ATT: return m_values[EST_ATT]*m_values[Q22_ATT];
    case Q23_REP: return m_values[EST_REP]*m_values[Q22_REP];
    case AST_ATT: return m_values[Q22_ATT]/m_values[Q11_ATT];
    case AST_REP: return m_values[Q22_REP]/m_values[Q11_REP];
    // Tabulated data
    default:      return m_values[type];
    }
}


void DebyeHuckleEvaluator::setDebyeLength(double Te, double ne)
{
    assert(Te >= 0.0);
    assert(ne >= 0.0);

    // Protect against small electron number densities
    ne = std::max(ne, 1.0e-16);

    // Including electron and ion contributions
    m_lambda = std::sqrt(0.5*EPS0*KB*Te/(ne*QE*QE));
}

void DebyeHuckleEvaluator::interpolate(double Tst)
{
    // Clip to table boundaries
    if (Tst <= sm_tstvec[0])
        m_values = sm_table.row(0);
    else if (Tst >= sm_tstvec[N_TST_POINTS-1])
        m_values = sm_table.row(N_TST_POINTS-1);
    else {
        // Find the index to interpolate between
        int i = 1; while (sm_tstvec[i] < Tst) ++i;

        // Interpolate
        m_values = (sm_table.row(i)-sm_table.row(i-1))*(Tst-sm_tstvec[i])/
            (sm_tstvec(i)-sm_tstvec(i-1)) + sm_table.row(i);
    }
}


// Initialize the T* index values
Array<double, N_TST_POINTS, 1> init_tstvec() {
    Array<double, N_TST_POINTS, 1> vec; vec <<
        0.1,   0.2,  0.3,  0.4,  0.6,  0.8,   1.0,   2.0,   3.0,  4.0,    6.0,   8.0,
       10.0,  20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0, 300.0, 400.0, 600.0, 800.0,
      1.0e3, 1.0e4;
    return vec;
}
Array<double, N_TST_POINTS, 1> DebyeHuckleEvaluator::sm_tstvec = init_tstvec();


// Initialize the table of collision integral values
Array<double, N_TST_POINTS, N_TABLE_TYPES, ColMajor> init_table() {
    Array<double, N_TST_POINTS, N_TABLE_TYPES, ColMajor> table; table <<
    // (T*)^2*Q11      (T*)^2*Q22      (T*)^2*Q14      (T*)^2*Q15      (T*)^2*Q24          B*              C*              E*
    //att     rep     att     rep     att     rep     att     rep     att     rep     att     rep     att     rep     att     rep
    0.0630, 0.0224, 0.0384, 0.0304, 0.0285, 0.0110, 0.0227, 0.0093, 0.0284, 0.0208, 1.4695, 1.3646, 0.7573, 0.7486, 0.8411, 0.8146,
    0.1364, 0.0511, 0.0967, 0.0697, 0.0460, 0.0221, 0.0353, 0.0181, 0.0652, 0.0445, 1.4577, 1.3865, 0.6585, 0.7104, 0.8121, 0.7836,
    0.1961, 0.0797, 0.1557, 0.1086, 0.0578, 0.0316, 0.0437, 0.0255, 0.0956, 0.0661, 1.4297, 1.3946, 0.6150, 0.6858, 0.7740, 0.7634,
    0.2480, 0.1072, 0.2100, 0.1459, 0.0669, 0.0399, 0.0500, 0.0317, 0.1208, 0.0856, 1.3987, 1.3976, 0.5852, 0.6675, 0.7460, 0.7483,
    0.3297, 0.1584, 0.3037, 0.2144, 0.0806, 0.0536, 0.0596, 0.0419, 0.1606, 0.1193, 1.3668, 1.3972, 0.5549, 0.6411, 0.7103, 0.7265,
    0.3962, 0.2050, 0.3818, 0.2757, 0.0910, 0.0648, 0.0668, 0.0500, 0.1914, 0.1476, 1.3425, 1.3933, 0.5356, 0.6221, 0.6878, 0.7108,
    0.4519, 0.2474, 0.4483, 0.3310, 0.0993, 0.0742, 0.0725, 0.0567, 0.2166, 0.1719, 1.3252, 1.3884, 0.5226, 0.6075, 0.6726, 0.6988,
    0.6467, 0.4177, 0.6840, 0.5460, 0.1269, 0.1065, 0.0915, 0.0792, 0.3007, 0.2587, 1.2798, 1.3627, 0.4901, 0.5632, 0.6344, 0.6628,
    0.7746, 0.5442, 0.8385, 0.6999, 0.1440, 0.1268, 0.1032, 0.0930, 0.3530, 0.3154, 1.2585, 1.3420, 0.4756, 0.5391, 0.6174, 0.6436,
    0.8719, 0.6455, 0.9541, 0.8197, 0.1566, 0.1416, 0.1118, 0.1030, 0.3914, 0.3574, 1.2442, 1.3255, 0.4663, 0.5230, 0.6073, 0.6311,
    1.0173, 0.8026, 1.1240, 1.0006, 0.1748, 0.1627, 0.1241, 0.1172, 0.4468, 0.4182, 1.2255, 1.3011, 0.4546, 0.5022, 0.5951, 0.6152,
    1.1259, 0.9230, 1.2486, 1.1385, 0.1880, 0.1777, 0.1330, 0.1272, 0.4869, 0.4620, 1.2133, 1.2833, 0.4471, 0.4886, 0.5878, 0.6052,
    1.2130, 1.0207, 1.3473, 1.2435, 0.1983, 0.1894, 0.1400, 0.1349, 0.5183, 0.4962, 1.2043, 1.2697, 0.4417, 0.4789, 0.5827, 0.5981,
    1.4972, 1.3431, 1.6626, 1.5892, 0.2307, 0.2254, 0.1616, 0.1589, 0.6163, 0.6027, 1.1798, 1.2295, 0.4273, 0.4529, 0.5696, 0.5797,
    1.6716, 1.5412, 1.8517, 1.7959, 0.2497, 0.2462, 0.1744, 0.1727, 0.6738, 0.6649, 1.1673, 1.2083, 0.4200, 0.4404, 0.5634, 0.5714,
    1.7984, 1.6847, 1.9872, 1.9438, 0.2634, 0.2610, 0.1835, 0.1825, 0.7149, 0.7089, 1.1590, 1.1945, 0.4154, 0.4326, 0.5595, 0.5663,
    1.9807, 1.8898, 2.1801, 2.1531, 0.2828, 0.2817, 0.1967, 0.1963, 0.7738, 0.7707, 1.1481, 1.1770, 0.4094, 0.4230, 0.5549, 0.5600,
    2.1123, 2.0368, 2.3184, 2.3019, 0.2969, 0.2964, 0.2062, 0.2061, 0.8164, 0.8145, 1.1409, 1.1659, 0.4056, 0.4170, 0.5521, 0.5561,
    2.2156, 2.1515, 2.4266, 2.4171, 0.3081, 0.3078, 0.2137, 0.2136, 0.8502, 0.8483, 1.1358, 1.1580, 0.4029, 0.4128, 0.5501, 0.5533,
    2.5427, 2.5087, 2.7672, 2.7713, 0.3433, 0.3429, 0.2373, 0.2370, 0.9566, 0.9534, 1.1220, 1.1365, 0.3956, 0.4013, 0.5448, 0.5456,
    2.7380, 2.7168, 2.9687, 2.9747, 0.3634, 0.3633, 0.2506, 0.2506, 1.0161, 1.0149, 1.1151, 1.1254, 0.3919, 0.3955, 0.5419, 0.5419,
    2.8780, 2.8635, 3.1121, 3.1177, 0.3778, 0.3778, 0.2602, 0.2602, 1.0588, 1.0583, 1.1106, 1.1181, 0.3894, 0.3920, 0.5401, 0.5398,
    3.0767, 3.0687, 3.3146, 3.3185, 0.3981, 0.3981, 0.2737, 0.2737, 1.1194, 1.1193, 1.1045, 1.1089, 0.3862, 0.3876, 0.5377, 0.5373,
    3.2185, 3.2135, 3.4583, 3.4610, 0.4125, 0.4125, 0.2833, 0.2833, 1.1628, 1.1624, 1.1005, 1.1033, 0.3841, 0.3849, 0.5361, 0.5358,
    3.3289, 3.3256, 3.5699, 3.5719, 0.4236, 0.4236, 0.2908, 0.2908, 1.1959, 1.1959, 1.0975, 1.0994, 0.3825, 0.3831, 0.5350, 0.5348,
    4.4759, 4.4763, 4.7211, 4.7211, 0.5388, 0.5388, 0.3675, 0.3675, 1.5413, 1.5413, 1.0734, 1.0733, 0.3702, 0.3702, 0.5265, 0.5265;
    return table;
}
Array<double, N_TST_POINTS, N_TABLE_TYPES, ColMajor>
    DebyeHuckleEvaluator::sm_table = init_table();

/**
 * Represents a collision integral computed with the Debye-Huckle potential
 * (Coulomb potential screened by the Debye length).
 */
class DebyeHuckleColInt : public CollisionIntegral
{
public:

    /**
     * Self registering constructor. Determines which integral to evaluate based
     * on the collision type.
     */
    DebyeHuckleColInt(CollisionIntegral::ARGS args) :
        CollisionIntegral(args)
    {
        assert(args.pair.type() == ATTRACTIVE || args.pair.type() == REPULSIVE);

        // Set the type of integral we need
        string kind = args.kind;
        CollisionType type = args.pair.type();

        if (kind == "Q11")
            m_type = (type == ATTRACTIVE ? Q11_ATT : Q11_REP);
        else if (kind == "Q12")
            m_type = (type == ATTRACTIVE ? Q12_ATT : Q12_REP);
        else if (kind == "Q13")
            m_type = (type == ATTRACTIVE ? Q13_ATT : Q13_REP);
        else if (kind == "Q14")
            m_type = (type == ATTRACTIVE ? Q14_ATT : Q14_REP);
        else if (kind == "Q15")
            m_type = (type == ATTRACTIVE ? Q15_ATT : Q15_REP);
        else if (kind == "Q22")
            m_type = (type == ATTRACTIVE ? Q22_ATT : Q22_REP);
        else if (kind == "Q23")
            m_type = (type == ATTRACTIVE ? Q23_ATT : Q23_REP);
        else if (kind == "Q24")
            m_type = (type == ATTRACTIVE ? Q24_ATT : Q24_REP);
        else if (kind == "Ast")
            m_type = (type == ATTRACTIVE ? AST_ATT : AST_REP);
        else if (kind == "Bst")
            m_type = (type == ATTRACTIVE ? BST_ATT : BST_REP);
        else if (kind == "Cst")
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

    double compute_(double T) { return sm_evaluator(T, m_type); }

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
ObjectProvider<DebyeHuckleColInt, CollisionIntegral> DebyeHuckle_ci("Debye-Huckel");

    } // namespace Transport
} // namespace Mutation
