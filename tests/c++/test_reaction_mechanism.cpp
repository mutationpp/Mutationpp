/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#include "mutation++.h"

#include <catch.hpp>
#include <Eigen/Dense>

#include <string>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Kinetics;
using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;
using namespace Catch;

/**
 * @brief Hard-coded reaction mechanism used for testing.
 *
 * The species are assumed to have the following order:
 *
 * e- N O NO NO+ N2 O2
 *
 * The reactions are:
 *
 * 1) N2 + M = 2N + M    kf = A1*sqrt(Th*Tv), kb = A1*Th/keq(Th)
 * 2) O2 + N = 2O + N    kf = A2*sqrt(Th*Tv), kb = A2*Th/keq(Th)
 * 3) O2 + e- = 2O + e-  kf = A3*Tv, kb = A3*Tv/keq(Tv)
 * 4) O2 + N = NO + O    kf = A4*Th, kb = A4*Th/keq(Th)
 * 5) N + O = NO+ + e-   kf = A5*Th, kb = A5*Tv/keq(Tv)
 */
class TestMechanism
{
public:

    /**
     * @brief Construct a new TestMechanism object.
     *
     * @param reversible Whether or not the reactions are reversible.
     */
    TestMechanism(bool reversible = true);

    /**
     * @brief Creates a mixture object representing this mechanism.
     *
     * @return SharedPtr<Mixture> A mixture object.
     */
    SharedPtr<Mixture> representativeMixture();

    /**
     * @brief Sets the current state of the mixture in the mechanism.
     *
     * @tparam Array Any type looking like an array.
     * @param densities Species densities (order = e- N O NO NO+ N2 O2).
     * @param Th Heavy translational and rotational temperature.
     * @param Tv Vibronic, electron translational temperature.
     */
    template <typename Array>
    void setState(Array densities, double Th, double Tv);

    /// Returns the forward rate coefficients.
    Eigen::ArrayXd& forwardRateCoefficients() { return m_kf; }
    /// Returns the reverse rate coefficients.
    Eigen::ArrayXd& reverseRateCoefficients() { return m_kb; }
    /// Returns the forward rates of progress.
    Eigen::ArrayXd& forwardRatesOfProgress() { return m_rf; }
    /// Returns the reverse rates of progress.
    Eigen::ArrayXd& reverseRatesOfProgress() { return m_rb; }
    /// Returns the reaction rates.
    Eigen::ArrayXd& netRatesOfProgress() { return m_rr; }
    /// Returns the production rates.
    Eigen::ArrayXd& productionRates() { return m_wdot; }
    /// Returns the production rate jacobians w.r.t. species densities.
    Eigen::MatrixXd& densityJacobian() { return m_jac_rho; }

private:

    bool m_reversible;
    std::string m_species;
    std::vector<std::string> m_formulas;
    Eigen::ArrayXd m_A;
    Eigen::ArrayXd m_Mw;
    Eigen::ArrayXd m_kf;
    Eigen::ArrayXd m_kb;
    Eigen::ArrayXd m_keq;
    Eigen::ArrayXd m_rf;
    Eigen::ArrayXd m_rb;
    Eigen::ArrayXd m_rr;
    Eigen::ArrayXd m_wdot;
    Eigen::MatrixXd m_jac_rho;
    SharedPtr<ThermoDB> m_thermo;

}; // class TestMechanism

TestMechanism::TestMechanism(bool reversible) :
    m_reversible(reversible), m_species("e- N O NO NO+ N2 O2"),
    m_A(5), m_Mw(7), m_kf(5), m_kb(5), m_keq(5), m_rf(5), m_rb(5), m_rr(5),
    m_wdot(7), m_jac_rho(7, 7)
{
    // Reaction formulas
    m_formulas.push_back("N2 + M = 2N + M");
    m_formulas.push_back("O2 + N = 2O + N");
    m_formulas.push_back("O2 + e- = 2O + e-");
    m_formulas.push_back("O2 + N = NO + O");
    m_formulas.push_back("N + O = NO+ + e-");

    if (!m_reversible) {
        std::vector<std::string>::iterator s;
        for (s = m_formulas.begin(); s != m_formulas.end(); ++s)
            s->insert(s->find('=')+1, ">");
    }

    // Rate constants
    for (int i = 0; i < 5; ++i)
        m_A(i) = (i+1) * 1000.0;

    // Species molecular weights
    double mwe = 0.00055;
    double mwn = 14.0067;
    double mwo = 15.9994;
    m_Mw << mwe, mwn, mwo, mwn + mwo, mwn + mwo - mwe, 2*mwn, 2*mwo;
    m_Mw /= 1000.0;

    // Thermodynamics database
    m_thermo = SharedPtr<ThermoDB>(Factory<ThermoDB>::create("RRHO", 0));
    m_thermo->load(m_species);
}

SharedPtr<Mixture> TestMechanism::representativeMixture()
{
    TemporaryFile mech_file(".xml");
    mech_file << "<mechanism>\n";
    for (int i = 0; i < m_formulas.size(); ++i) {
        mech_file
            << "<reaction formula=\"" << m_formulas[i] <<  "\" >\n"
            << "    <arrhenius A=\"" << m_A(i) << "\" n=\"1\" T=\"0.0\" />\n"
            << "</reaction>\n";
    }
    mech_file  << "</mechanism>\n";
    mech_file.close();

    MixtureOptions mixture_options;
    mixture_options.setStateModel("ChemNonEqTTv");
    mixture_options.setThermodynamicDatabase("RRHO");
    mixture_options.setSpeciesDescriptor(m_species);
    mixture_options.setDefaultComposition("N2: 0.79, O2: 0.21");
    mixture_options.setMechanism(mech_file.filename());

    return SharedPtr<Mixture>(new Mixture(mixture_options));
}

template <typename Array>
void TestMechanism::setState(Array densities, double Th, double Tv)
{
    // e- N O NO NO+ N2 O2
    // 0  1 2 3  4   5  6

    m_kf(0) = m_A(0) * std::sqrt(Th * Tv);
    m_kf(1) = m_A(1) * std::sqrt(Th * Tv);
    m_kf(2) = m_A(2) * Tv;
    m_kf(3) = m_A(3) * Th;
    m_kf(4) = m_A(4) * Th;

    m_thermo->gibbs(
        Th, Th, Th, Th, Th, m_thermo->standardPressure(),
        m_wdot.data(), NULL, NULL, NULL, NULL);
    m_keq(0) = ONEATM / (RU * Th) * std::exp(-2.0*m_wdot(1) + m_wdot(5));
    m_keq(1) = ONEATM / (RU * Th) * std::exp(-2.0*m_wdot(2) + m_wdot(6));
    m_keq(3) = std::exp(-m_wdot(3) - m_wdot(2) + m_wdot(6) + m_wdot(1));

    m_thermo->gibbs(
        Tv, Tv, Tv, Tv, Tv, m_thermo->standardPressure(),
        m_wdot.data(), NULL, NULL, NULL, NULL);
    m_keq(2) = ONEATM / (RU * Tv) * std::exp(-2.0*m_wdot(2) + m_wdot(6));
    m_keq(4) = std::exp(-m_wdot(4) - m_wdot(0) + m_wdot(1) + m_wdot(2));

    m_kb.setZero();
    if (m_reversible) {
        m_kb(0) = m_A(0) * Th / m_keq(0);
        m_kb(1) = m_A(1) * Th / m_keq(1);
        m_kb(2) = m_A(2) * Tv / m_keq(2);
        m_kb(3) = m_A(3) * Th / m_keq(3);
        m_kb(4) = m_A(4) * Tv / m_keq(4);
    }

    for (int i = 0; i < 7; ++i)
        m_wdot[i] = densities[i] / m_Mw[i];
    double M = m_wdot.tail(6).sum();

    m_rf(0) = m_kf(0) * m_wdot(5) * M;
    m_rf(1) = m_kf(1) * m_wdot(6) * m_wdot(1);
    m_rf(2) = m_kf(2) * m_wdot(6) * m_wdot(0);
    m_rf(3) = m_kf(3) * m_wdot(6) * m_wdot(1);
    m_rf(4) = m_kf(4) * m_wdot(1) * m_wdot(2);
    m_rr = m_rf;

    m_rb(0) = m_kb(0) * m_wdot(1) * m_wdot(1) * M;
    m_rb(1) = m_kb(1) * m_wdot(2) * m_wdot(2) * m_wdot(1);
    m_rb(2) = m_kb(2) * m_wdot(2) * m_wdot(2) * m_wdot(0);
    m_rb(3) = m_kb(3) * m_wdot(3) * m_wdot(2);
    m_rb(4) = m_kb(4) * m_wdot(4) * m_wdot(0);
    if (m_reversible) m_rr -= m_rb;

    // Compute jacobian terms
    Eigen::VectorXd drr(7), stoich(7);
    m_jac_rho.setConstant(0.0);

    // 1) N2 + M = 2N + M
    drr.setConstant(0.0);
    drr(1) = -2*m_kb(0)*m_wdot(1); // N
    drr(5) = m_kf(0); // N2
    stoich << 0, 2, 0, 0, 0, -1, 0;
    m_jac_rho += stoich * (M * drr).transpose();
    for (int j = 1; j < 7; ++j)
        m_jac_rho.col(j) += stoich * m_rr(0) / M;

    // 2) O2 + N = 2O + N
    drr.setConstant(0.0);
    drr(1) = m_kf(1)*m_wdot(6) - m_kb(1)*m_wdot(2)*m_wdot(2); // N
    drr(2) = -2*m_kb(1)*m_wdot(1)*m_wdot(2); // O
    drr(6) = m_kf(1)*m_wdot(1); // O2
    stoich << 0, 0, 2, 0, 0, 0, -1;
    m_jac_rho += stoich * drr.transpose();

    // 3) O2 + e- = 2O + e-
    drr.setConstant(0.0);
    drr(0) = m_kf(2)*m_wdot(6) - m_kb(2)*m_wdot(2)*m_wdot(2);
    drr(2) = -2*m_kb(2)*m_wdot(2)*m_wdot(0);
    drr(6) = m_kf(2)*m_wdot(0);
    stoich << 0, 0, 2, 0, 0, 0, -1;
    m_jac_rho += stoich * drr.transpose();

    // 4) O2 + N = NO + O
    drr.setConstant(0.0);
    drr(1) = m_kf(3)*m_wdot(6);
    drr(2) = -m_kb(3)*m_wdot(3);
    drr(3) = -m_kb(3)*m_wdot(2);
    drr(6) = m_kf(3)*m_wdot(1);
    stoich << 0, -1, 1, 1, 0, 0, -1;
    m_jac_rho += stoich * drr.transpose();

    // 5) N + O = NO+ + e-
    drr.setConstant(0.0);
    drr(0) = -m_kb(4)*m_wdot(4);
    drr(1) = m_kf(4)*m_wdot(2);
    drr(2) = m_kf(4)*m_wdot(1);
    drr(4) = -m_kb(4)*m_wdot(0);
    stoich << 1, -1, -1, 0, 1, 0, 0;
    m_jac_rho += stoich * drr.transpose();

    m_wdot(0) = m_rr(4);
    m_wdot(1) = 2*m_rr(0) - m_rr(3) - m_rr(4);
    m_wdot(2) = 2*m_rr(1) + 2*m_rr(2) + m_rr(3) - m_rr(4);
    m_wdot(3) = m_rr(3);
    m_wdot(4) = m_rr(4);
    m_wdot(5) = -m_rr(0);
    m_wdot(6) = -m_rr(1) - m_rr(2) - m_rr(3);

    m_wdot *= m_Mw;
    
    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j)
            m_jac_rho(i,j) *= (m_Mw(i)/m_Mw(j));
}

/// Jumps to another permutation of x.
bool stillPermuting(std::vector<double>& x)
{
    bool permuting = true;
    for (int i = 0; i < 100; ++i)
        permuting &= std::next_permutation(x.begin(), x.end());
    return permuting;
}

void testMech(bool isRev)
{
    TestMechanism mech(isRev);
    SharedPtr<Mixture> mix = mech.representativeMixture();

    // Create a vector of "densities" in sorted order
    std::vector<double> rho(mix->nSpecies()+1);
    for (int i = 0; i < rho.size(); ++i)
        rho[i] = i*0.1;

    std::vector<double> T(2, 1000.0);
    Eigen::ArrayXd test_nr(mix->nReactions());
    Eigen::ArrayXd test_ns(mix->nSpecies());
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> 
        test_jac(mix->nSpecies(), mix->nSpecies());

    // Loop over all permutations of the densities and temperatures
    do {
        for (int i = 0; i < 10; ++i) {
            T[0] = 500.0 + i*1000.0;
            for (int j = 0; j < 10; ++j) {
                T[1] = 500.0 + j*1000.0;

                mix->setState(rho.data(), T.data(), 1);
                mech.setState(rho.data(), T[0], T[1]);

                mix->forwardRateCoefficients(test_nr.data());
                REQUIRE(test_nr.isApprox(mech.forwardRateCoefficients()));

                mix->backwardRateCoefficients(test_nr.data());
                REQUIRE(test_nr.isApprox(mech.reverseRateCoefficients()));

                mix->forwardRatesOfProgress(test_nr.data());
                REQUIRE(test_nr.isApprox(mech.forwardRatesOfProgress()));

                mix->backwardRatesOfProgress(test_nr.data());
                REQUIRE(test_nr.isApprox(mech.reverseRatesOfProgress()));

                mix->netRatesOfProgress(test_nr.data());
                REQUIRE(test_nr.isApprox(mech.netRatesOfProgress()));

                mix->netProductionRates(test_ns.data());
                REQUIRE(test_ns.isApprox(mech.productionRates()));

                mix->jacobianRho(test_jac.data());
                REQUIRE(test_jac.isApprox(mech.densityJacobian()));
            }
        }
    } while (stillPermuting(rho));
}

/**
 * Tests if the TestMechanism evaluation matches Kinetics.
 */
TEST_CASE("Test mechanism is evaluated correctly", "[kinetics]")
{
    testMech(true);
    testMech(false);
}
