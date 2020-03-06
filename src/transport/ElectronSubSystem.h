/**
 * @file ElectronsubSystem.h
 *
 * Provides functions for computing the electron transport properties.
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

#ifndef ELECTRON_SUB_SYSTEM_H
#define ELECTRON_SUB_SYSTEM_H

#include "Thermodynamics.h"
#include "CollisionDB.h"
#include "Constants.h"

#include <Eigen/Dense>
#include <vector>

namespace Mutation {
    namespace Transport {

/**
 * Provides functions which solve the electron-heavy transport systems.  It is
 * convenient to package all of these functions into a separate class.
 */
class ElectronSubSystem
{
public:
    /// Constructor.
    ElectronSubSystem(
        Mutation::Thermodynamics::Thermodynamics& thermo,
        CollisionDB& collisions) :
            m_thermo(thermo),
            m_collisions(collisions),
            m_alpha(thermo.nHeavy()),
            m_alpha_B(thermo.nHeavy(), 3),
            m_chi2(thermo.nHeavy()),
            m_chi2_B(thermo.nHeavy(), 3)
    { }

    /// Destructor.
    ~ElectronSubSystem() { }

    /// Isotropic electric conductivity in S/m.
    double electricConductivity(int order = 3);

    /// Anisotropic electric conductivity in S/m.
    Eigen::Vector3d electricConductivityB(int order = 3);

    /// Isotropic electron diffusion coefficient.
    double electronDiffusionCoefficient(int order = 3);

    /// Anisotropic electron diffusion coefficient.
    Eigen::Vector3d electronDiffusionCoefficientB(int order = 3);

    /// Isotropic second-order electron diffusion coefficient.
    double electronDiffusionCoefficient2(
        const Eigen::Ref<const Eigen::MatrixXd>& Dij, int order = 3);

    /// Anisotropic second-order electron diffusion coefficient.
    Eigen::Vector3d electronDiffusionCoefficient2B(
        const Eigen::Ref<const Eigen::MatrixXd>& Dij, int order = 3);

    /// Isotropic electron thermal conductivity in W/m-K.
    double electronThermalConductivity(int order = 3);

    /// Anisotropic electron thermal conductivity in W/m-K.
    Eigen::Vector3d electronThermalConductivityB(int order = 3);

    /// Isotropic alpha coefficients.
    const Eigen::VectorXd& alpha(int order = 3);

    /// Anisotropic alpha coefficients.
    const Eigen::Matrix<double,-1,3>& alphaB(int order = 3);

    /// Isotropic electron thermal diffusion ratio.
    double electronThermalDiffusionRatio(int order = 3);

    /// Anisotropic electron thermal diffusion ratio.
    Eigen::Vector3d electronThermalDiffusionRatioB(int order = 3);

    /// Isotropic second-order electron thermal diffusion ratios.
    const Eigen::VectorXd& electronThermalDiffusionRatios2(int order = 3);

    /// Anisotropic second-order electron thermal diffusion ratios.
    const Eigen::Matrix<double,-1,3>& electronThermalDiffusionRatios2B(int order = 3);


    /**
     * Returns the factor which must be multiplied by the result of the
     * Lee() function to retrieve \f$\L_{ee}^{pq}\f$.
     */
    double Leefac() {
        const double Te = m_thermo.Te();
        const double me = m_collisions.mass()[0];
        return 16.*m_thermo.P()/(3.*KB*Te)*std::sqrt(me/(TWOPI*KB*Te));
    }

    /**
     * Returns the symmetric \f$\Lambda_{ee}\f$ matrix of the given size divided
     * by the \f$\frac{64x_e}{75k_B}\sqrt{\frac{m_e}{2\pi k_B T_e}}\f$.
     */
    template <int SIZE>
    Eigen::Matrix<double, SIZE, SIZE> Lee();

    /**
     * Returns the \f$\L_{ee}^{Bpq}\f$ matrix with the given size.
     */
    template <int SIZE>
    Eigen::Matrix<double, SIZE, SIZE> LBee();


protected:

    /**
     * Helper function to evaluate \f$ sum_{i\in\mathcal{H}} x_i a_i \f$.
     */
    template <typename E>
    double dotxh(const Eigen::ArrayBase<E>& a) {
        return (m_collisions.X()*a).tail(m_collisions.nHeavy()).sum();
    }

    template <int P>
    Eigen::Matrix<double,P,P> L1inv() {
        return Lee<P>().inverse()/Leefac();
    }

//    template <int P>
//    Eigen::Matrix<std::complex<double>,P,P> L2inv()
//    {
//        Eigen::Matrix<std::complex<double>,P,P> L2;
//        L2.real() = Leefac()*Lee<P>();
//        L2.imag() = LBee<P>();
//        return L2.inverse();
//    }

    template <int P>
    Eigen::Vector3d electronThermalConductivityB();

    template <int P>
    double electronDiffusionCoefficient() {
        return L1inv<P>()(0,0);
    }

    template <int P>
    Eigen::Vector3d electronDiffusionCoefficientB()
    {
        // Linear system matrices
        Eigen::Matrix<double,P,P> L1 = Leefac() * Lee<P>();
        Eigen::Matrix<std::complex<double>,P,P> L2;
        L2.real() = L1;
        L2.imag() = LBee<P>();

        Eigen::Vector3d De;
        De(0) = L1.inverse().eval()(0);

        std::complex<double> alpha = L2.inverse().eval()(0);
        De(1) = alpha.real();
        De(2) = alpha.imag();

        return De;
    }

    template <int P>
    const Eigen::VectorXd& alpha();

    template <int P>
    const Eigen::Matrix<double,-1,3>& alphaB();

    template <int P>
    double electronThermalDiffusionRatio();

    template <int P>
    Eigen::Vector3d electronThermalDiffusionRatioB();

    template <int P>
    const Eigen::VectorXd& electronThermalDiffusionRatios2();

    template <int P>
    const Eigen::Matrix<double,-1,3>& electronThermalDiffusionRatios2B();

//    template <int P, typename RHS>
//    Eigen::Matrix<double, P, 1> solveRealSysP0(const RHS& rhs) {
//        return Lee<XI>().inverse() * rhs<XI>();
//    }
//
//    template <int XI, typename RHS>
//    Eigen::Matrix<std::complex<double>, XI, 1> solveImagSysP0(const RHS& rhs)
//    {
//        Eigen::Matrix<std::complex<double>, XI, XI> sys;
//        sys.real(Lee<XI>()); sys.imag(LBee<XI>());
//        return sys.inverse() * rhs<XI>();
//    }

private:

    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    CollisionDB m_collisions;

    Eigen::VectorXd m_alpha;
    Eigen::Matrix<double,-1,3> m_alpha_B;
    Eigen::VectorXd m_chi2;
    Eigen::Matrix<double,-1,3> m_chi2_B;

}; // class ElectronSubSystem


template <int SIZE>
Eigen::Matrix<double, SIZE, SIZE> ElectronSubSystem::Lee()
{
    Eigen::Matrix<double, SIZE, SIZE> L;
    const int nh = m_collisions.nHeavy();
    const double xe = m_collisions.X()(0);

    // L00ee
    const Eigen::ArrayXd& Q11 = m_collisions.Q11ei();
    L(0,0) = dotxh(Q11);

    // Return if done
    if (SIZE == 1) return L;

    const Eigen::ArrayXd& Q12 = m_collisions.Q12ei();
    const Eigen::ArrayXd& Q13 = m_collisions.Q13ei();
    const double Q22 = m_collisions.Q22ee();

    // L01ee, L10ee
    L(0,1) = dotxh(2.5*Q11-3.*Q12);
    L(1,0) = L(0,1);

    // L11ee
    L(1,1) = dotxh(6.25*Q11-15.*Q12+12.*Q13);
    L(1,1) += xe*SQRT2*Q22;

    // Return if done
    if (SIZE == 2) return L;

    const Eigen::ArrayXd& Q14 = m_collisions.Q14ei();
    const Eigen::ArrayXd& Q15 = m_collisions.Q15ei();
    const double Q23 = m_collisions.Q23ee();
    const double Q24 = m_collisions.Q24ee();

    // L02ee, L20ee
    L(0,2) = dotxh(35./8.*Q11-10.5*Q12+6.*Q13);
    L(2,0) = L(0,2);

    // L12ee, L21ee
    L(1,2) = dotxh(175./16.*Q11-315./8.*Q12+57.*Q13-30.*Q14);
    L(1,2) += xe*SQRT2*(1.75*Q22 - 2.*Q23);
    L(2,1) = L(1,2);

    // L22ee
    L(2,2) = dotxh(1225./64.*Q11-735./8.*Q12+199.5*Q13-210.*Q14+90.*Q15);
    L(2,2) += xe*SQRT2*(77./16.*Q22-7.*Q23+5.*Q24);

    return L;
}

template <int SIZE>
Eigen::Matrix<double, SIZE, SIZE> ElectronSubSystem::LBee()
{
    Eigen::Matrix<double, SIZE, SIZE> LB;
    const double fac = m_thermo.getBField()*QE/(KB*m_thermo.Te());

    LB(0,0) = fac;
    if (SIZE == 1) return LB;

    LB(1,1) = 2.5*fac;
    if (SIZE == 2) return LB;

    LB(2,2) = 35./8.*fac;
    return LB;
}

template <int P>
Eigen::Vector3d ElectronSubSystem::electronThermalConductivityB()
{
    const double fac = 2.5*m_thermo.numberDensity()*m_thermo.X()[0]*KB;

    Eigen::Matrix<double,P-1,P-1> L1 =
        Leefac()*Lee<P>().template bottomRightCorner<P-1,P-1>();
    Eigen::Matrix<std::complex<double>,P-1,P-1> L2;
    L2.real() = L1;
    L2.imag() = LBee<P>().template bottomRightCorner<P-1,P-1>();

    Eigen::Vector3d lambda;
    lambda(0) = L1.inverse().eval()(0);

    std::complex<double> alpha = L2.inverse().eval()(0);
    lambda(1) = alpha.real();
    lambda(2) = alpha.imag();

    return 2.5*fac*lambda;
}


/// Small helper class which provides the \f$\beta^{pD_i}\f$.
template <int P>
class BetaDi
{
public:
    BetaDi(
        const Mutation::Thermodynamics::Thermodynamics& thermo,
        CollisionDB& collisions) : m_beta(P, thermo.nHeavy())
    {
        const double fac = 16.0/3.0*thermo.numberDensity()*
            std::sqrt(collisions.mass()(0)/(TWOPI*KB*thermo.Te()));
        //const double fac = 16.0/3.0*thermo.numberDensity()/std::sqrt(TWOPI*thermo.Te());
        const int nh = thermo.nHeavy();

        const Eigen::ArrayXd& Q11 = collisions.Q11ei();
        m_beta.row(0) = fac*(collisions.X()*Q11).tail(nh);

        if (P == 1) return;

        const Eigen::ArrayXd& Q12 = collisions.Q12ei();
        m_beta.row(1) = fac*(collisions.X()*(2.5*Q11 - 3.0*Q12)).tail(nh);

        if (P == 2) return;

        const Eigen::ArrayXd& Q13 = collisions.Q13ei();
        m_beta.row(2) =
            fac*(collisions.X()*(4.375*Q11 - 10.5*Q12 + 6.0*Q13)).tail(nh);
    }

    const Eigen::Matrix<double,P,Eigen::Dynamic>& operator() () {
        return m_beta;
    }

private:

    Eigen::Matrix<double,P,Eigen::Dynamic> m_beta;

};

template <int P>
const Eigen::VectorXd& ElectronSubSystem::alpha()
{
    // Compute the system solution
    BetaDi<P> beta(m_thermo, m_collisions);
    Eigen::Matrix<double,P,P> L1inv = Lee<P>().inverse();

    for (int i = 0; i < m_thermo.nHeavy(); ++i)
        m_alpha(i) = (L1inv.col(0)).dot(beta().col(i));

    return (m_alpha /= Leefac());
}

template <int P>
const Eigen::Matrix<double,-1,3>& ElectronSubSystem::alphaB()
{
    // Linear system matrices
    Eigen::Matrix<double,P,P> L1 = Leefac() * Lee<P>();
    Eigen::Matrix<std::complex<double>,P,P> L2;
    L2.real() = L1;
    L2.imag() = LBee<P>();

    // Inverses
    Eigen::Matrix<double,P,P> L1inv = L1.inverse();
    Eigen::Matrix<std::complex<double>,P,P> L2inv = L2.inverse();

    // Compute the system solution
    BetaDi<P> beta(m_thermo, m_collisions);
    std::complex<double> sol;
    Eigen::Matrix<double,P,1> b;

    for (int i = 0; i < m_thermo.nHeavy(); ++i) {
        b = beta().col(i);
        m_alpha_B(i,0) = (L1inv.col(0)).dot(b);

        sol = (L2inv.col(0)).dot(b);
        m_alpha_B(i,1) = sol.real();
        m_alpha_B(i,2) = sol.imag();
    }

    return m_alpha_B;
}

template <int P>
double ElectronSubSystem::electronThermalDiffusionRatio()
{
    Eigen::Matrix<double,P,P> L1 = Lee<P>();

    return 2.5 * (L1.template topRightCorner<1,P-1>() *
        L1.template bottomRightCorner<P-1,P-1>().inverse().template leftCols<1>())(0);
}

template <int P>
Eigen::Vector3d ElectronSubSystem::electronThermalDiffusionRatioB()
{
    // Linear system matrices
    Eigen::Matrix<double,P,P> L1 = Leefac() * Lee<P>();
    Eigen::Matrix<std::complex<double>,P,P> L2;
    L2.real() = L1;
    L2.imag() = LBee<P>();

    Eigen::Vector3d chi;
    chi(0) = (L1.template topRightCorner<1,P-1>() *
        L1.template bottomRightCorner<P-1,P-1>().inverse().template leftCols<1>())(0);

    std::complex<double> sol = (L2.template topRightCorner<1,P-1>() *
        L2.template bottomRightCorner<P-1,P-1>().inverse().template leftCols<1>())(0);
    chi(1) = sol.real();
    chi(2) = sol.imag();

    return 2.5*chi;
}

template <int P>
const Eigen::VectorXd& ElectronSubSystem::electronThermalDiffusionRatios2()
{
    BetaDi<P> beta(m_thermo, m_collisions);
    m_chi2 = -2.5 * (L1inv<P>() * beta()).row(1);
    return m_chi2;
}

template <int P>
const Eigen::Matrix<double,-1,3>& ElectronSubSystem::electronThermalDiffusionRatios2B()
{
    // Linear system matrices
    Eigen::Matrix<double,P,P> L1 = Leefac() * Lee<P>();
    Eigen::Matrix<std::complex<double>,P,P> L2;
    L2.real() = L1;
    L2.imag() = LBee<P>();

    BetaDi<P> beta(m_thermo, m_collisions);
    m_chi2_B.col(0) = -2.5 * (L1.inverse() * beta()).row(1);

    Eigen::VectorXcd sol = (L2.inverse() * beta()).row(1);
    m_chi2_B.col(1) = -2.5 * sol.real();
    m_chi2_B.col(2) = -2.5 * sol.imag();

    return m_chi2_B;
}


    } // namespace Transport
} // namespace Mutation

#endif // ELECTRON_SUB_SYSTEM_H
