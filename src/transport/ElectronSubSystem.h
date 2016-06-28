/**
 * @file ElectronsubSystem.h
 *
 * Provides functions for computing the electron transport properties.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2016 von Karman Institute for Fluid Dynamics (VKI)
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
            m_alpha(thermo.nHeavy())
    { }

    /// Destructor.
    ~ElectronSubSystem() { }

    /// Isotropic electric conductivity in S/m.
    double electricConductivity(int order = 2);

    /// Anisotropic electric conductivity in S/m.
    Eigen::Vector3d electricConductivityB(int order = 2);

    /// Isotropic electron thermal conductivity in W/m-K.
    double electronThermalConductivity(int order = 2);

    /// Anisotropic electron thermal conductivity in W/m-K.
    Eigen::Vector3d electronThermalConductivityB(int order = 2);

    const Eigen::VectorXd& alpha(int order = 2);

protected:

    /**
     * Returns the factor which must be multiplied by the result of the
     * Lee() function to retrieve \f$\L_{ee}^{pq}\f$.
     */
    double Leefac() {
        const double Te = m_thermo.Te();
        const double me = m_collisions.mass()[0];
        const double xe = m_thermo.X()[0];
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

    /**
     * Helper function to evaluate \f$ sum_{i\in\mathcal{H}} x_i a_i \f$.
     */
    template <typename E>
    double dotxh(const Eigen::ArrayBase<E>& a) {
        return (m_collisions.X()*a).tail(m_collisions.nHeavy()).sum();
    }

    template <int XI>
    Eigen::Vector3d sigmaB()
    {
        const double fac =
            m_thermo.numberDensity()*m_thermo.X()[0]*QE*QE/(KB*m_thermo.Te());

        Eigen::Vector3d sigma;
        sigma(0) = m_collisions.L1inv<XI>()(0,0);

        std::complex<double> alpha = m_collisions.L2inv<XI>()(0,0);
        sigma(1) = alpha.real();
        sigma(2) = alpha.imag();

        return fac*sigma;
    }

    template <int ORDER>
    Eigen::Matrix<double,ORDER,ORDER> L1inv() {
        return Lee<ORDER>().inverse()/Leefac();
    }

    template <int ORDER>
    Eigen::Matrix<std::complex<double>,ORDER,ORDER> L2inv()
    {
        Eigen::Matrix<std::complex<double>,ORDER,ORDER> L2 = LBee<ORDER>();
        L2.real() += Leefac()*Lee<ORDER>();
        return L2.inverse();
    }

    template <int P>
    Eigen::Vector3d electronThermalConductivityB()
    {
        const double fac = 2.5*m_thermo.numberDensity()*m_thermo.X()[0]*KB;

        Eigen::Matrix<double,P-1,P-1> L1 =
            m_collisions.Leefac()*
            m_collisions.Lee<P>().template block<P-1,P-1>(1,1);
        Eigen::Matrix<double,P-1,P-1> L1inv = L1.inverse();
        Eigen::Matrix<std::complex<double>,P-1,P-1> L2 =
            m_collisions.LBee<P>().template block<P-1,P-1>(1,1);
        L2.real() += L1;
        Eigen::Matrix<std::complex<double>,P-1,P-1> L2inv = L2.inverse();

        Eigen::Vector3d lambda;
        lambda(0) = 2.5*L1inv(0,0);

        std::complex<double> alpha = 2.5*L2inv(0,0);
        lambda(1) = alpha.real();
        lambda(2) = alpha.imag();

        return fac*lambda;
    }

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
    const double Q22 = m_collisions.Q22ee()(0);

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


    } // namespace Transport
} // namespace Mutation

#endif // ELECTRON_SUB_SYSTEM_H
