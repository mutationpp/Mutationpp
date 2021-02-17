/**
 * @file Transport.cpp
 *
 * @brief Implements Transport class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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


#include "Constants.h"
#include "DiffusionMatrix.h"
#include "ThermalConductivityAlgorithm.h"
#include "Transport.h"
#include "ViscosityAlgorithm.h"

#include <iostream>
#include <Eigen/Dense>

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::Config;
using namespace Eigen;

namespace Mutation {
    namespace Transport {

using Mutation::Thermodynamics::Thermodynamics;

//==============================================================================

Transport::Transport(
    Thermodynamics& thermo, const std::string& viscosity, const std::string& lambda)
    : m_thermo(thermo),
      m_collisions("collisions.xml", thermo),
      mp_esubsyst(NULL),
      mp_viscosity(NULL),
      mp_thermal_conductivity(NULL),
      mp_diffusion_matrix(NULL),
      mp_wrk1(NULL),
      mp_tag(NULL)
{
    // Setup the electron subsystem object
    mp_esubsyst = new ElectronSubSystem(m_thermo, m_collisions);

    // Load the viscosity calculator
    setViscosityAlgo(viscosity);

    // Load the thermal conductivity calculator
    setThermalConductivityAlgo(lambda);

    // Load the diffusion matrix calculator
    setDiffusionMatrixAlgo("Ramshaw");

    // Allocate work array storage
    mp_wrk1 = new double [m_thermo.nGas()*3];
    mp_wrk2 = mp_wrk1 + m_thermo.nGas();
    mp_wrk3 = mp_wrk2 + m_thermo.nGas();
    if(m_thermo.nEnergyEqns() > 1) {
        mp_tag  = new int [m_thermo.nEnergyEqns()*5];
        m_thermo.getTagModes(mp_tag);
    }

    //thermo.stateModel()->notifyOnUpdate(this);
}

//==============================================================================

Transport::~Transport()
{
    delete mp_esubsyst;
    delete mp_viscosity;
    delete mp_thermal_conductivity;
    delete mp_diffusion_matrix;

    delete [] mp_wrk1;
    delete [] mp_tag;
}

//==============================================================================

void Transport::setViscosityAlgo(const std::string& algo)
{
    if (mp_viscosity != NULL)
        delete mp_viscosity;

    try {
        mp_viscosity = Factory<ViscosityAlgorithm>::create(algo, m_collisions);
    } catch (Error& e) {
        e << "\nWas trying to set the viscosity algorithm.";
        throw;
    }
}

//==============================================================================

double Transport::viscosity() { return mp_viscosity->viscosity(); }

//==============================================================================

void Transport::setThermalConductivityAlgo(const std::string& algo)
{
    if (mp_thermal_conductivity != NULL)
        delete mp_thermal_conductivity;

    try {
        mp_thermal_conductivity =
            Factory<ThermalConductivityAlgorithm>::create(algo, m_collisions);
    } catch (Error& e) {
        e << "\nWas trying to set the thermal conductivity algorithm.";
        throw;
    }
}

//==============================================================================

double Transport::heavyThermalConductivity()
{
    return mp_thermal_conductivity->thermalConductivity();
}

//==============================================================================

void Transport::setDiffusionMatrixAlgo(const std::string& algo)
{
    if (mp_diffusion_matrix != NULL)
        delete mp_diffusion_matrix;

    try {
        mp_diffusion_matrix =
            Factory<DiffusionMatrix>::create(algo, m_collisions);
    } catch (Error& e) {
        e << "\nWas trying to set the diffusion matrix algorithm.";
        throw;
    }
}

//==============================================================================

const Eigen::MatrixXd& Transport::diffusionMatrix()
{
    return mp_diffusion_matrix->diffusionMatrix();
}

//==============================================================================

void Transport::heavyThermalDiffusionRatios(double* const p_k)
{
    mp_thermal_conductivity->thermalDiffusionRatios(p_k);
}

//==============================================================================

void Transport::frozenThermalConductivityVector(double* const p_lambda)
{
	const int neq = m_thermo.nEnergyEqns();
    double lambda_th, lambda_te, lambda_rot, lambda_vib, lambda_elec;

    if(neq > 1) {
        lambda_th   = heavyThermalConductivity();
        lambda_te   = electronThermalConductivity();
        lambda_rot  = rotationalThermalConductivity();
        lambda_vib  = vibrationalThermalConductivity();
        lambda_elec = electronicThermalConductivity();

        for (int i_eq = 0; i_eq < neq; ++i_eq) {
           p_lambda[i_eq] = mp_tag[i_eq*5+0] * lambda_th
                          + mp_tag[i_eq*5+1] * lambda_te
                          + mp_tag[i_eq*5+2] * lambda_rot
                          + mp_tag[i_eq*5+3] * lambda_vib
                          + mp_tag[i_eq*5+4] * lambda_elec;
        }
    } else {
        p_lambda[0] = frozenThermalConductivity();
    }

}

//==============================================================================

double Transport::internalThermalConductivity(double T)
{
    m_thermo.thermoDB()->cpint(T, mp_wrk1);
    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nGas()));
}


//==============================================================================

double Transport::rotationalThermalConductivity()
{
	m_thermo.speciesCpOverR(
        m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(), m_thermo.Tel(),
        NULL, NULL, mp_wrk1, NULL, NULL);

    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nGas()));
}

//==============================================================================

double Transport::vibrationalThermalConductivity()
{
	m_thermo.speciesCpOverR(
        m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(), m_thermo.Tel(),
        NULL, NULL, NULL, mp_wrk1, NULL);

    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nGas()));
}

//==============================================================================

double Transport::electronicThermalConductivity()
{
	m_thermo.speciesCpOverR(
        m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(), m_thermo.Tel(),
        NULL, NULL, NULL, NULL, mp_wrk1);

    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nGas()));
}

//==============================================================================

double Transport::reactiveThermalConductivity()
{
	// Compute dX_i/dT
    m_thermo.dXidT(mp_wrk1);

    // Compute the thermal diffusion ratios
    heavyThermalDiffusionRatios(mp_wrk2);

    // @todo Add the electron thermal diffusion ratios

    // Combine to get the driving forces
    for (int i = 0; i < m_thermo.nGas(); i++)
        mp_wrk1[i] += mp_wrk2[i] / m_thermo.T();

    // Compute the diffusion velocities
    double E;
    stefanMaxwell(mp_wrk1, mp_wrk2, E);

    // Unitless enthalpies
    m_thermo.speciesHOverRT(mp_wrk1);

    double lambda = 0.0;
    //const double* const X = m_thermo.X();
    const double rho = m_thermo.density();
    const double* const Y = m_thermo.Y();
    for (int i = 0; i < m_thermo.nGas(); ++i)
        lambda -= mp_wrk1[i] / m_thermo.speciesMw(i) * mp_wrk2[i] * Y[i] * rho;

    return (RU * m_thermo.T() * lambda);
}

//==============================================================================

double Transport::butlerBrokawThermalConductivity()
{
    const int ns = m_thermo.nGas();
    const int ne = m_thermo.nElements();
    const int nr = ns - ne;
    const int nh = m_thermo.nHeavy();
    const int a  = ns - nh;

    const Eigen::MatrixXd& E = m_thermo.elementMatrix();

    // Find the elements
    Eigen::VectorXi ei(ne);
    for (int i = 0; i < ne; ++i)
        ei(i) = m_thermo.speciesIndex(m_thermo.elementName(i));

    // Generate formation reactions (assuming elements are always included)
    Eigen::MatrixXd nu(nr, ns);
    nu.setConstant(0);
    for (int i = 0, ir = 0; i < ns; ++i) {
        if (E.row(i).array().abs().sum() == 1)
            continue;

        // Add formation reaction for the species
        for (int j = 0; j < ne; ++j)
            nu(ir,ei(j)) = -E(i,j);
        nu(ir++,i) = 1;
    }

    // Compute the Delta H vector
    Eigen::VectorXd H(ns);
    Eigen::VectorXd DH(nr);
    m_thermo.speciesHOverRT(H.data());
    DH = (nu * H) * KB * m_thermo.T();

    // Compute the system matrix
    Eigen::MatrixXd A(nr,nr); A.setConstant(0.0);
    Eigen::MatrixXd nDij(ns,ns);
    if (m_thermo.hasElectrons()) {
        for (int i = 0; i < ns; ++i) {
            nDij(i,0) = m_collisions.nDei()(i);
            nDij(0,i) = m_collisions.nDei()(i);
        }
    }
    for (int i = a, s = 0; i < ns; ++i) {
        for (int j = i; j < ns; ++j, s++) {
            nDij(i,j) = m_collisions.nDij()(s);
            nDij(j,i) = m_collisions.nDij()(s);
        }
    }
    const Eigen::ArrayXd X = m_collisions.X().max(1.0e-16);

    // Heavy species
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j <= i; ++j) {
            for (int k = 0; k < ns - 1; ++k) {
                for (int l = k + 1; l < ns; ++l) {
                    A(i,j) += X(k)*X(l)/nDij(k,l)*
                              (nu(i,k)/X(k)-nu(i,l)/X(l))*
                              (nu(j,k)/X(k)-nu(j,l)/X(l));
                }
            }
        }
    }

    //A *= RU*m_thermo.T()*m_thermo.numberDensity()/m_thermo.P();

    Eigen::VectorXd W = Eigen::LDLT<Eigen::MatrixXd, Eigen::Lower>(A).solve(DH);

    return W.dot(DH)/(KB*m_thermo.T()*m_thermo.T());
}

//==============================================================================

double Transport::soretThermalConductivity()
{
	// @todo This is super inefficient, should fix
	Eigen::VectorXd work(m_thermo.nGas());

	// Compute dX_i/dT
    m_thermo.dXidT(mp_wrk1);

    // Compute the thermal diffusion ratios
    heavyThermalDiffusionRatios(mp_wrk2);

    // @todo Add the electron thermal diffusion ratios

    // Combine to get the driving forces
    for (int i = 0; i < m_thermo.nGas(); i++)
        mp_wrk1[i] += mp_wrk2[i] / m_thermo.T();

    // Compute the diffusion velocities
    double E;
    stefanMaxwell(mp_wrk1, work.data(), E);

    double lambda = 0.0;
    for (int i = 0; i < m_thermo.nGas(); i++)
        lambda -= mp_wrk2[i]*work[i];

    return (m_thermo.P()*lambda);
}

//==============================================================================

void Transport::equilDiffFluxFacs(double* const p_F)
{
	// Get some state data
	const int ns = m_thermo.nGas();
	const int ne = m_thermo.nElements();
	const double* const p_Y = m_thermo.Y();
	const double* const p_X = m_thermo.X();
	const double rho = m_thermo.density();
	const double p   = m_thermo.P();

	const MatrixXd& Dij = diffusionMatrix();
	const Eigen::MatrixXd& nu  = m_thermo.elementMatrix();

	for (int i = 0; i < ns; ++i) {
		mp_wrk2[i] = 0.0;
		for (int j = 0; j < ns; ++j)
			mp_wrk2[i] += Dij(i,j)*mp_wrk1[j];
		mp_wrk2[i] *= -rho*p_Y[i]/m_thermo.speciesMw(i);
	}

	for (int k = 0; k < ne; ++k) {
		p_F[k] = 0.0;
		double mwk = m_thermo.atomicMass(k);
		for (int i = 0; i < ns; ++i)
			p_F[k] += nu(i,k)*mwk*mp_wrk2[i];
	}

	m_thermo.speciesHOverRT(mp_wrk1);
	p_F[ne] = 0.0;
	for (int i = 0; i < ns; ++i)
		p_F[ne] += mp_wrk1[i]*mp_wrk2[i];
	p_F[ne] *= RU * m_thermo.T();
}

//==============================================================================

void Transport::equilDiffFluxFacsP(double* const p_F)
{
	const int ns = m_thermo.nGas();
	const double p = m_thermo.P();
	const double* const p_Y = m_thermo.Y();
	const double* const p_X = m_thermo.X();

    // Compute the dXj/dP term
    m_thermo.dXidP(mp_wrk1);

	// Add the (x_j - y_j)/p term
	for (int i = 0; i < ns; ++i)
		mp_wrk1[i] += (p_X[i] - p_Y[i])/p;

    // Compute the element averaged diffusion coefficients
	equilDiffFluxFacs(p_F);
}

//==============================================================================

void Transport::equilDiffFluxFacsT(double* const p_F)
{
	const int ns = m_thermo.nGas();
	const double T  = m_thermo.T();

    Eigen::Map<Eigen::ArrayXd> work1(mp_wrk1, ns);
    Eigen::Map<Eigen::ArrayXd> work2(mp_wrk2, ns);

    // Compute the dXj/dT term
    m_thermo.dXidT(work1.data());

	// Add thermal diffusion ratio term
	heavyThermalDiffusionRatios(work2.data());
    work1 += work2 / T;

    // @todo Add electron thermal diffusion ratios

    // Compute the element averaged diffusion coefficients
	equilDiffFluxFacs(p_F);
}

//==============================================================================

void Transport::equilDiffFluxFacsZ(double* const p_F)
{
	const int ns = m_thermo.nGas();
	const int ne = m_thermo.nElements();
	const double* const p_X = m_thermo.X();
    const double T   = m_thermo.T();
    const double p   = m_thermo.P();

	// Loop over each element
	for (int l = 0; l < ne; ++l) {
	   // Compute the dXj/dZl term using a finite difference
	   m_thermo.elementFractions(p_X, mp_wrk1);
	   double h = std::max(mp_wrk1[l]*1.0e-6, 1.0e-10);
	   mp_wrk1[l] += h;
	   m_thermo.equilibriumComposition(T, p, mp_wrk1, mp_wrk2);

	   for (int i = 0; i < ns; ++i)
		   mp_wrk1[i] = (mp_wrk2[i]-p_X[i])/h;

	   // Compute the element averaged diffusion coefficients
	   equilDiffFluxFacs(p_F + l*(ne+1));
	}

	// Be sure to set the state back in the equilibrium solver in case other
	// calculations rely on the correct element potential values
	m_thermo.elementFractions(p_X, mp_wrk1);
	m_thermo.equilibriumComposition(T, p, mp_wrk1, mp_wrk2);
}

//==============================================================================

void Transport::stefanMaxwell(
    const double* const p_dp, double* const p_V, double& E, int order)
{
    stefanMaxwell(m_thermo.T(), m_thermo.Te(), p_dp, p_V, E, order);
}


void Transport::stefanMaxwell(double Th, double Te,
    const double* const p_dp, double* const p_V, double& E, int order)
{
    const int ns = m_thermo.nGas();
    const int k  = ns - m_thermo.nHeavy();
    const double nd = m_thermo.numberDensity();

    static ArrayXd X;
    X = Map<const ArrayXd>(m_thermo.X(), ns) + 1.0e-16;
    X /= X.sum();

    ArrayXd qi(ns), Mi(ns);
    for (int i = 0; i < ns; ++i) {
        qi(i) = m_thermo.speciesCharge(i);
        Mi(i) = m_thermo.speciesMw(i);
    }

    // Form the SM matrix
    MatrixXd G = MatrixXd::Zero(ns+k,ns+k);

    // electron subsystem
    double s = 0.0;
    double a = 0.0;
    if (k == 1) {
        const ArrayXd& nDei = m_collisions.nDei();
        ArrayXd phi; smCorrectionsElectron(order, phi);
        for (int i = 1; i < ns; ++i) {
            const double fac = Te/Th*X(0)*X(i)/nDei(i)*nd*(1.0+phi(i));
            G(0,0) += fac;
            G(i,0) = -fac;
            G(i,i) = Te/Th*fac;
            G(0,i) = -fac;
        }
        G(0,0) *= Th/Te;

        ArrayXd kappa = X*(qi - Mi * (qi*X).sum() / (Mi*X).sum()) / (KB*Th);

        s = kappa.matrix().norm();
        G.row(ns).head(ns) = kappa / s;
        G.col(ns).head(ns) = kappa / s;
        G(0,ns) *= Th/Te;

        a = nDei.mean();
    }

    // heavy subsystem
    const ArrayXd& nDij = m_collisions.nDij();
    ArrayXd phi; smCorrectionsHeavy(order, phi);
    for (int i = k, is = 1; i < ns; ++i, ++is) {
        for (int j = i+1; j < ns; ++j, ++is) {
            const double fac = X(i)*X(j)/nDij(is)*nd*(1.0+phi(is));
            G(i,j) = -fac;
            G(i,i) += fac;
            G(j,i) = -fac;
            G(j,j) += fac;
        }
    }

    // add mass constraint
    //a = nd / std::max(a, nDij.maxCoeff());
    a = G.diagonal().maxCoeff();
    ArrayXd Y = X * Mi / (Mi*X).sum();
    G.topLeftCorner(ns,ns) += a * Y.matrix() * Y.matrix().transpose();

    // Form the right hand side
    VectorXd b = VectorXd::Zero(ns+k);
    b.head(ns) = -Map<const VectorXd>(p_dp, ns);
    if (k == 1) b(0) *= Th/Te;

    // Solve the system
    // Modified solver because old one failing
    VectorXd x = G.colPivHouseholderQr().solve(b);
   // VectorXd x = G.householderQr().solve(b);

    // Retrieve the solution
    Map<VectorXd>(p_V, ns) = x.head(ns);
    E = (k == 1 ? x(ns)/s : 0.0);

    // Apply Ramshaw to fix roundoff errors
    Map<ArrayXd>(p_V, ns) -= (Map<ArrayXd>(p_V, ns)*Map<const ArrayXd>(m_thermo.Y(), ns)).sum();
}

void Transport::smCorrectionsElectron(int order, Eigen::ArrayXd& phi)
{
    const int ns = m_thermo.nGas();
    phi = ArrayXd::Zero(ns);

    if (order == 1)
        return;

    const ArrayXd X = Map<const ArrayXd>(m_thermo.X(), ns).max(1.0e-16);
    const ArrayXd& nDei = m_collisions.nDei();
    const Matrix2d Lee = mp_esubsyst->Lee<2>();
    const ArrayXd& L01 = m_collisions.L01ei();

    phi = 25./4.*KB*nDei/(X*X(0))*Lee(0,1)/Lee(1,1)*L01;
    std::cout << phi << "\n" << std::endl;
}

void Transport::smCorrectionsHeavy(int order, Eigen::ArrayXd& phi)
{
    //order = 2;
    //std::cout << "smCorrectionsHeavy: order = " << order << std::endl;

    const int nh = m_thermo.nHeavy();
    phi = ArrayXd::Zero(nh*(nh+1)/2);

    if (order == 1)
        return;

    // Use second order corrections
    const int ns = m_thermo.nGas();
    const int k  = ns - nh;
    ArrayXd X = Map<const ArrayXd>(m_thermo.X(), ns) + 1.0e-16;
    X /= X.sum();

    const ArrayXd& mi = m_collisions.mass();
    const ArrayXd& Ast = m_collisions.Astij();
    const ArrayXd& Bst = m_collisions.Bstij();
    const ArrayXd& Cst = m_collisions.Cstij();
    const ArrayXd& nDij = m_collisions.nDij();
    const ArrayXd& etai = m_collisions.etai();

    // Compute the Lam01 matrix
    MatrixXd Lam01(nh,nh);
    double fac;
    Lam01.diagonal().setZero();
    for (int j = 0, si = 1; j < nh; ++j, ++si) {
        for (int i = j+1; i < nh; ++i, ++si) {
            fac = X(i+k)*X(j+k)/(mi(i+k)+mi(j+k))*(12.*Cst(si)-10.)/(25.*KB*nDij(si));
            Lam01(i,j) = fac*mi(i+k);
            Lam01(j,i) = fac*mi(j+k);
            Lam01(i,i) -= Lam01(j,i);
            Lam01(j,j) -= Lam01(i,j);
        }
    }


    // Compute the Glamh matrix (only lower part)
    MatrixXd Glamh(nh,nh);
    Glamh.diagonal().array() = 4./(15.0*KB)*(X*X*mi).tail(nh)/etai;
    double miij, mjij;
    int ik, jk;
    for (int j = 0, si = 1; j < nh; ++j, ++si) {
        jk = j+k;
        for (int i = j+1; i < nh; ++i, ++si) {
            ik = i+k;
            miij = mi(ik) / (mi(ik) + mi(jk));
            mjij = mi(jk) / (mi(ik) + mi(jk));
            fac  = X(ik)*X(jk) / (nDij(si)*25.*KB);
            Glamh(i,j) = fac * miij * mjij *
                (16.0 * Ast(si) + 12.0 * Bst(si) - 55.0);
            Glamh(i,i) += fac * (miij * (30.0 * miij + 16.0 * mjij *
                Ast(si)) + mjij * mjij * (25.0 - 12.0 * Bst(si)));
            Glamh(j,j) += fac * (mjij * (30.0 * mjij + 16.0 * miij *
                Ast(si)) + miij * miij * (25.0 - 12.0 * Bst(si)));
        }
    }

    //std::cout << Glamh << std::endl << std::endl;
    //std::cout << Lam01 << std::endl << std::endl;

    // Get LDLT factorization of the Glamh matrix
    LDLT<MatrixXd, Lower> ldlt(Glamh);

    // Compute second order corrections
    VectorXd beta(nh);
    for (int j = 0, is = 1; j < nh; ++j, ++is) {
        beta = ldlt.solve(Lam01.row(j).transpose());
        for (int i = j+1; i < nh; ++i, ++is)
            phi(is) = nDij(is)/(X(i+k)*X(j+k)) * beta.dot(Lam01.row(i));
    }
    phi *= 6.25 * KB;

    //std::cout << phi << "\n" << std::endl;
}


//void Transport::stefanMaxwell(
//    const double* const p_dp, double* const p_V, double& E)
//{
//	const int ns = m_thermo.nGas();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//
//    // Electron offset
//    const int k = m_thermo.hasElectrons() ? 1 : 0;
//    const int nh = ns-k;
//
//    Map<const ArrayXd> dp(p_dp, ns);
//    Map<ArrayXd> V(p_V, ns), qi(mp_wrk3, ns);
//
////    // Need to place a tolerance on X and Y
//    const double tol = 1.0e-16;
//    static ArrayXd X, Y; Y.resize(ns); // only gets resized the first time
//    X = Map<const ArrayXd>(m_thermo.X(), ns).max(tol); // Place a tolerance on X
//    m_thermo.convert<X_TO_Y>(&X[0], &Y[0]);
////    Eigen::Map<const Eigen::ArrayXd> X = m_collisions.X();
////    Eigen::Map<const Eigen::ArrayXd> Y = m_collisions.Y();
//
//    // Get reference to binary diffusion coefficients
//    const ArrayXd& nDij = m_collisions.nDij();
//
//    // Compute mixture charge
//    for (int i = 0; i < ns; ++i)
//        qi[i] = m_thermo.speciesCharge(i);
//    const double q = (X*qi).sum();
//
//    // Compute kappas (store in V)
//    V = (X*qi - q*Y)/(KB*Th);
//
//    // Compute ambipolar electric field
//    E = (k > 0 ? dp[0]/V[0] : 0.0);
//
//    // Compute right hand side
//    static VectorXd b; b.resize(nh);
//    b.array() = -dp.tail(nh) + E*V.tail(nh);
//
//    // Compute singular system matrix (lower triangular part)
//    double fac;
//    static MatrixXd G; G.resize(nh,nh);
//    G.triangularView<Lower>() = MatrixXd::Zero(nh, nh);
//
//    for (int j = 0, si = 1; j < nh; ++j, ++si) {
//        for (int i = j+1; i < nh; ++i, ++si) {
//            fac = X[i+k]*X[j+k]/nDij(si)*nd;
//            G(i,i) += fac;
//            G(j,j) += fac;
//            G(i,j) = -fac;
//        }
//    }
//
//    // Compute the constraints that are applied to the matrix
//    fac = X[0]*qi[0];
//    V.tail(nh) = Y.tail(nh);
//    if (k > 0) V.tail(nh) -= X.tail(nh)*qi.tail(nh)*Y[0]/fac;
//
//    // Add mass balance relation to make matrix nonsigular
//    G.selfadjointView<Lower>().rankUpdate(
//        V.matrix().tail(nh), nd/nDij.mean());
//
//    static Eigen::LDLT<MatrixXd, Lower> ldlt;
//    ldlt.compute(G);
//    V.tail(nh) = ldlt.solve(b).array();
//
//    // Compute electron diffusion velocity if need be
//    if (k == 0)
//        return;
//
//    V[0] = -(X.tail(nh)*qi.tail(nh)*V.tail(nh)).sum() / fac;
//}

//==============================================================================


double Transport::meanFreePath()
{
    // Thermo properties
         const int ns = m_thermo.nGas();
         const double Th = m_thermo.T();
           const double Te = m_thermo.Te();
                 const double nd = m_thermo.numberDensity();
                    const double* const X = m_thermo.X();
                    const double me = m_thermo.speciesMw(0)/NA;
                   const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
        const double Q11ee = m_collisions.Q11ee();
        const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();


                                         double sum = 0.0;
            sum +=X[0]*X[0]*Q11ee;
                for (int i = 1; i < ns; ++i)
            {
            sum +=X[i]*X[0]*Q11ei(i);
            ;}

                                           for (int i = 1; i < ns; ++i)
                                                for (int j = 1; j < ns; ++j)
                                                               sum += X[i]*X[j]*Q11(i,j);

                                                                   return 1.0/(nd*sum);
                                                                 }
//==============================================================================


double Transport::electronMeanFreePath()
{
    if (!m_thermo.hasElectrons())
        return 0.0;

    const int ns = m_thermo.nGas();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();

    const double Q11ee = m_collisions.Q11ee();
    const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
    double sum = 0.0;
            sum +=X[0]*X[0]*Q11ee;
                        for (int i = 1; i < ns; ++i)
                        {
                        sum +=X[i]*X[0]*Q11ei(i);
                        ;}


    return 1.0/(nd*sum);
}

//==============================================================================
double Transport::speciesThermalSpeed(const int& i) const
{
	if (i < m_thermo.hasElectrons()){
	    const double T = m_thermo.Te();
        return sqrt(8.0*RU*T/(PI*m_thermo.speciesMw(i)));
	}

    const double T = m_thermo.T();
    return sqrt(8.0*RU*T/(PI*m_thermo.speciesMw(i)));
}

//==============================================================================
double Transport::averageHeavyThermalSpeed()
{
    const int ns = m_thermo.nGas();
    const int nh = m_thermo.nHeavy();
    const double Th = m_thermo.T();
    const double* const X = m_thermo.X();

    double ave_mw = 0.0;
    for (int i = ns-nh; i < ns; ++i)
        ave_mw += X[i] * m_thermo.speciesMw(i);

    return sqrt(8.0*RU*Th/(PI*ave_mw));

}

//==============================================================================
double Transport::electronThermalSpeed()
{
    const double Te = m_thermo.Te();
    const double Me = m_thermo.speciesMw(0);
    return sqrt(8*RU*Te/(PI*Me));
}

//==============================================================================
double Transport::averageHeavyCollisionFreq()
{
  return meanFreePath()/averageHeavyThermalSpeed();
}
//==============================================================================
double Transport::electronHeavyCollisionFreq()
{
    if (!m_thermo.hasElectrons())
        return 0.0;
    return electronThermalSpeed()/electronMeanFreePath();
}
//==============================================================================
double Transport::coulombMeanCollisionTime()
{
    if (!m_thermo.hasElectrons())
        return 0.0;

    const int ns = m_thermo.nGas();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    double sum = 0.0;
    const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
    const double Q11ee = m_collisions.Q11ee();
    const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
                        sum +=X[0]*X[0]*Q11ee;
                        for (int i = 1; i < ns; ++i)
                        {
                        sum +=X[i]*X[0]*Q11ei(i);
                        ;}

                                           for (int i = 1; i < ns; ++i)
                                                for (int j = 1; j < ns; ++j)
                                                               sum += X[i]*X[j]*Q11(i,j);

     return (3.0/16.0)*1.0/(nd*sum)/averageHeavyThermalSpeed();//electronThermalSpeed();
     }

//==============================================================================
double Transport::hallParameter()
{
    if (!m_thermo.hasElectrons())
        return 0.0;

    const double me = m_thermo.speciesMw(0)/NA;
    const double B = m_thermo.getBField();

    return QE*B*electronMeanFreePath()/(me*electronThermalSpeed());
}
//==============================================================================

// double Transport::parallelDiffusionCoefficient()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;
// //first order
// //return fac2/lam00;

// //second order
// //return fac2/(lam00-lam01*lam01/lam11);
// // third order
// return -fac2*(lam12*lam12-lam11*lam22)/(lam00*lam11*lam22-lam00*lam12*lam12-lam22*lam01*lam01+2.0*lam01*lam02*lam12-lam11*lam02*lam02);
// }
// //==============================================================================
// double Transport::perpDiffusionCoefficient()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25./4.*nd*KB*(Th/(X[0]*Te)+(1.-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;
// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;
// lamB00 = QE*B/(KB*fac3*Te); //fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;

// //first order
// //return fac2*lam00/(lam00*lam00+lamB00*lamB00);


// //second order

// // double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
//   // double numerator = lam11*(lam00*lam11-lamB00*lamB11-lam01*lam01) + lamB11*(lamB00*lam11+lamB11*lam00);
//     // return fac2*numerator/denominator;



// //third order
// double denominator = (lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02);
//     double numerator = -((lam12*lam12-lam11*lam22+lamB11*lamB22)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) - (lam11*lamB22+lam22*lamB11)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02));


// return fac2*(numerator/denominator);

// }
// //==============================================================================
// double Transport::transverseDiffusionCoefficient()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;
// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;


// // First order
// //     //return -fac2*lamB00/(lam00*lam00+lamB00*lamB00);


// // Second order
// // double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
//  //           double numerator = lamB11*(lam00*lam11-lamB00*lamB11-lam01*lam01) - lam11*(lamB00*lam11+lamB11*lam00);
// //                    return fac2*numerator/denominator;


// //third order

// double denominator = (lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02);
//     double numerator = ((lam11*lamB22+lam22*lamB11)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam12*lam12-lam11*lam22+lamB11*lamB22)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02));
// return fac2*(numerator/denominator);
// }
// //==============================================================================
// double Transport::parallelThermalDiffusionCoefficient()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;
// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;


// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;

// // First order
//  // return fac2*(5.0/2.0)/lam00;

//    //Second order
// //        return fac2*(5.0/2.0)*lam01/(lam00*lam11-lam01*lam01);
// // Third order
// return fac2*(-5.0/2.0)*(-lam12/(lam11*lam22-lam12*lam12));
// }

// //==============================================================================
// double Transport::perpThermalDiffusionCoefficient()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;

// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;
// //first order
// //return fac2*(5.0/2.0)*lam00/(lam00*lam00+lamB00*lamB00);

// //second order
// //double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
// //  double numerator = -lam01*(lam00*lam11-lamB00*lamB11-lam01*lam01);

// //third order
// double numerator = lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11);
// double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);
// //
// //        //return fac2*(-5.0/2.0)*numerator/denominator;

// return fac2*(-5.0/2.0)*(numerator/denominator);
// }

// //==============================================================================
// double Transport::transverseThermalDiffusionCoefficient()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;

// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;
// //first order
// // return -fac2*(-5.0/2.0)*lamB00/(lam00*lam00+lamB00*lamB00);


// //second order
// //double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
// //    double numerator = lam01*(lamB00*lam11+lamB11*lam00);
// //return fac2*(-5.0/2.0)*(numerator/denominator);

// // Third order

// double numerator = lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11);
// double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);

//              return fac2*(-5.0/2.0)*numerator/denominator;
// }
// //==============================================================================
// double Transport::sigmaParallel()
// {
//     if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//         return 0.0;
// //avant il y avait pas le fac2
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double* const X = m_thermo.X();
//     const double nd = m_thermo.numberDensity();
//     const double P = m_thermo.P();
//     double fac2=(4.*X[0]/(25.*nd*KB));

//     return fac2*(nd*KB*Te/(P))*(4.0/25.0)*(X[0]*QE)*(X[0]*QE)*(1.0/(KB*KB*Te))*(parallelDiffusionCoefficient()/(fac2));

// }
// //==============================================================================
// double Transport::sigmaPerpendicular()
// {
//     if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//         return 0.0;
// //avant il y avait pas le fac2
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double* const X = m_thermo.X();
//     const double nd = m_thermo.numberDensity();
//     const double P = m_thermo.P();
//     double fac2=(4.*X[0]/(25.*nd*KB));

//     return fac2*(nd*KB*Te/(P))*(4.0/25.0)*(X[0]*QE)*(X[0]*QE)*(1.0/(KB*KB*Te))*(perpDiffusionCoefficient()/(fac2));

// }
// //==============================================================================
// double Transport::sigmaTransverse()
// {
//     if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//         return 0.0;
// //avant il y avait pas le fac2
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double* const X = m_thermo.X();
//     const double nd = m_thermo.numberDensity();
//     const double P = m_thermo.P();
//     double fac2=(4.*X[0]/(25.*nd*KB));

//     return fac2*(nd*KB*Te/(P))*(4.0/25.0)*(X[0]*QE)*(X[0]*QE)*(1.0/(KB*KB*Te))*(transverseDiffusionCoefficient()/(fac2));

// }

// //==============================================================================
// //
// double Transport::parallelElectronThermalConductivity()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;


// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;

// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;


// //second order
// //return (X[0]*X[0])/(lam11);
// //third order
// return (X[0]*X[0]*lam22)/(lam11*lam22-lam12*lam12);

// }


// //==============================================================================
// double Transport::perpElectronThermalConductivity()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;

// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;
//  // First order (pas sur que ca marche)
//  //     //return (X[0]*X[0])*lam00/(lam00*lam00+lamB00*lamB00);
// //second order
// //return (X[0]*X[0])*lam11/((lam11)*(lam11)+(lamB11)*(lamB11));
// //third order
// //
// double numerator = lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11);
// double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);
// ////
// return (X[0]*X[0])*numerator/denominator;
// //
// //
// //

// }

// //==============================================================================
// double Transport::transverseElectronThermalConductivity()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=L(0,0)/fac;
// lam01=L(0,1)/fac;
// lam02=L(0,2)/fac;
// lam11=L(1,1)/fac;
// lam12=L(1,2)/fac;
// lam22=L(2,2)/fac;


// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;

// //second order
// //return -(X[0]*X[0])*lamB11/((lam11)*(lam11)+(lamB11)*(lamB11));
// //third order
// double numerator = lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11);
// double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);
// //
//  return (X[0]*X[0])*numerator/denominator;





// }

// //==============================================================================
// double Transport::parallelThermalDiffusionRatio()
// {

// return (parallelThermalDiffusionCoefficient())/(parallelDiffusionCoefficient());}

// //==============================================================================
// double Transport::transverseThermalDiffusionRatio()
// {

// return (perpDiffusionCoefficient()*transverseThermalDiffusionCoefficient()-transverseDiffusionCoefficient()*perpThermalDiffusionCoefficient())/(perpDiffusionCoefficient()*perpDiffusionCoefficient()+transverseDiffusionCoefficient()*transverseDiffusionCoefficient());
// }
// //==============================================================================
// double Transport::perpThermalDiffusionRatio()
// {

// return (perpDiffusionCoefficient()*perpThermalDiffusionCoefficient()+transverseThermalDiffusionCoefficient()*transverseDiffusionCoefficient())/(perpDiffusionCoefficient()*perpDiffusionCoefficient()+transverseDiffusionCoefficient()*transverseDiffusionCoefficient());}

// //==============================================================================
// double Transport::taueLambda()
// {
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
// const double Th = m_thermo.T();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const int ns = m_thermo.nGas();
// const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
// const double Q11ee = m_collisions.Q11ee();
// const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
// const Eigen::ArrayXd& Q12ei = m_collisions.Q12ei();
// const Eigen::ArrayXd& Q13ei = m_collisions.Q13ei();
// const double Q22ee = m_collisions.Q22ee();
// double tauelambda=0.0;
// double Ve0=0.0;
// double collisionalsum=0.0;
// Ve0=std::sqrt(KB*Te/(2*PI*me));
// collisionalsum += nd*X[0]*(SQRT2*Q22ee);
// collisionalsum += nd*X[1]*(25.0/4.0*Q11ei(1)-15.0*Q12ei(1)+12.0*Q13ei(1));
// tauelambda=1.0/((64.0/75.0)*Ve0*(collisionalsum));

// return tauelambda;
// }
// //==============================================================================
// //////====== THis transport braginskii coefficient is only for a fully ionized plasma like e-, A+ with an order one for the ionization
// //

// double Transport::taueLambdaBr()
// {
//  double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const double B = m_thermo.getBField();
// const double Th = m_thermo.T();

// double lam = 0.0;
// double taue = 0.0;
// double temporary = 0.0;
// double TT=0.0;
// double nne=0.0;
// double mme=0.0;
// double omegae= 0.0;

// TT=Te*8.62E-5;
// nne=nd*X[0]*1.0E-6;
// temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
// if (TT>36.2){
// temporary=25.25+std::log(TT)-0.5*std::log(nne);
// }
// taue=3.1616*3.44E5*TT*std::sqrt(TT)/(temporary*nne);

// return taue;
// }

// //==============================================================================
// ////////====== THis is the coefficient friction from mutation ++
// ////
// double Transport::coefficientFriction()
// {
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
// const double Th = m_thermo.T();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const int ns = m_thermo.nGas();
// const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
// const double Q11ee = m_collisions.Q11ee();
// const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
// const Eigen::ArrayXd& Q12ei = m_collisions.Q12ei();
// const Eigen::ArrayXd& Q13ei = m_collisions.Q13ei();
// const double Q22ee = m_collisions.Q22ee();
// double frictioncoeff=0.0;
// double A=0.0;
// double BB=0.0;
// A=(5.0/2.0)*Q11ei(1)-3.0*Q12ei(1);
// BB=25.0/4.0*Q11ei(1)-15.0*Q12ei(1)+12.0*Q13ei(1)+SQRT2*Q22ee;
// frictioncoeff=(5.0/2.0)*(Q11ei(1)*A-(A)*(A)*(A)/BB)/(Q11ei(1)*BB-A*A);
// // BE careful the friction coefficient is only for a fully ionized mixture composed of hydrogen+ and electron , mixture Hp !!!!!
// //
// return frictioncoeff;
// }

// //==============================================================================
// double Transport::perpfriccoeffBr()
// {
//  double nd  = m_thermo.numberDensity();
//   const double Te = m_thermo.Te();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double* const X = m_thermo.X();
//    const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();

//    double perpfriccoeffBr=0.0;
//    double lam = 0.0;
//    double taue = 0.0;
//    double temporary = 0.0;
//    double TT=0.0;
//    double nne=0.0;
//    double mme=0.0;
//    double omegae= 0.0;
//    double beta0=0.0;
//    double xx=0.0;
//    double delta=0.0;
//    omegae=QE*B/me;
//    TT=Te*8.62E-5;
//    nne=nd*X[0]*1.0E-6;
//    temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//    if (TT>36.2){
//    temporary=25.25+std::log(TT)-0.5*std::log(nne);
//    }
//    taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//    xx=omegae*taue;
//    delta=xx*xx*xx*xx+14.8*xx*xx+3.77;

//    perpfriccoeffBr=(5.1*xx*xx+2.7)/delta;
//    return perpfriccoeffBr;
//   }
// //==============================================================================
// double Transport::transfriccoeffBr()
// {
//  double nd  = m_thermo.numberDensity();
//   const double Te = m_thermo.Te();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double* const X = m_thermo.X();
//    const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();

//    double transfriccoeffBr=0.0;
//    double lam = 0.0;
//    double taue = 0.0;
//    double temporary = 0.0;
//    double TT=0.0;
//    double nne=0.0;
//    double mme=0.0;
//    double omegae= 0.0;
//    double beta0=0.0;
//    double xx=0.0;
//    double delta=0.0;
//    omegae=QE*B/me;
//    TT=Te*8.62E-5;
//    nne=nd*X[0]*1.0E-6;
//    temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//    if (TT>36.2){
//    temporary=25.25+std::log(TT)-0.5*std::log(nne);
//    }
//    taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//    xx=omegae*taue;
//    delta=xx*xx*xx*xx+14.8*xx*xx+3.77;

//    transfriccoeffBr=xx*(1.5*xx*xx+3.053)/delta;
//    return transfriccoeffBr;
//   }




// //==============================================================================

// double Transport::tauViscosity()
// {
// const double mh = m_thermo.speciesMw(1)/NA;
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
// const double Th = m_thermo.T();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const int ns = m_thermo.nGas();
// const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
// const Eigen::ArrayXd& Q22 = m_collisions.Q22ij();
// const double Q11ee = m_collisions.Q11ee();
// const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
// const Eigen::ArrayXd& Q12ei = m_collisions.Q12ei();
// const Eigen::ArrayXd& Q13ei = m_collisions.Q13ei();
// const double Q22ee = m_collisions.Q22ee();

// double tauvisc=0.0;
// double Vh0=0.0;
// double collisionsum=0.0;

// Vh0=std::sqrt(KB*Th/(PI*mh));
// collisionsum=nd*X[1]*(4.8*Q22(0)+(5.0/3.0)*Q11(0));

// tauvisc=1/(Vh0*collisionsum);

//  return tauvisc;
// }

// //==============================================================================
// double Transport::tauViscosityBr()
// {
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const double B = m_thermo.getBField();
// const double Th = m_thermo.T();
// double taui = 0.0;
//   double temporary = 0.0;
//   double TT=0.0;
//   double nne=0.0;
// TT=Th*8.62E-5;
//   nne=nd*X[0]*1.0E-6;
//   temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
// if (TT>36.2){
// temporary=25.25+std::log(TT)-0.5*std::log(nne);
// }

//   taui=0.96*2.09E7*TT*std::sqrt(TT)/(nne*temporary);
//   return taui;

// }

// //==============================================================================
// double Transport::tauLambdaHeavy()
// {
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//   const double mh = m_thermo.speciesMw(1)/NA;
//   const double* const X = m_thermo.X();
//   const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();
// const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
// const Eigen::ArrayXd& Q22 = m_collisions.Q22ij();
// const double Q11ee = m_collisions.Q11ee();
// const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
// const Eigen::ArrayXd& Q12ei = m_collisions.Q12ei();
// const Eigen::ArrayXd& Q13ei = m_collisions.Q13ei();
// const Eigen::ArrayXd& Bstar = m_collisions.Bstij();
// const double Q22ee = m_collisions.Q22ee();
//    double taulambdah=0.0;
//    double Vh0=0.0;
//    double collisionsum=0.0;
//  Vh0=std::sqrt(KB*Th/(PI*mh));
//   collisionsum=nd*X[1]*((16.0/600.0)*Q11ei(1)*(55.0-12.0*Bstar(1)+16.0*(Q22(1,1))/(Q11(1,1)))+(64.0/75.0)*Q22(1,1));
//   taulambdah=1/(Vh0*collisionsum);

//   return taulambdah;
// }
// //==============================================================================
// //
// double Transport::tauLambdaHeavyBr()
// {
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const double B = m_thermo.getBField();
// const double Th = m_thermo.T();
// double taui = 0.0;
//   double temporary = 0.0;
//   double TT=0.0;
//   double nne=0.0;
// TT=Th*8.62E-5;
//     nne=nd*X[0]*1.0E-6;
//     temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
// if (TT>36.2){
// temporary=25.25+std::log(TT)-0.5*std::log(nne);
// }
//     taui=3.906*2.09E7*TT*std::sqrt(TT)/(nne*temporary);
//     return taui;

// }

// //==============================================================================
// double Transport::tauEnergy()
// {
// double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//   const double mh = m_thermo.speciesMw(1)/NA;
// const double me =m_thermo.speciesMw(0)/NA;
//   const double* const X = m_thermo.X();
//   const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();
// const Eigen::ArrayXd& Q11 = m_collisions.Q11ij();
// const Eigen::ArrayXd& Q22 = m_collisions.Q22ij();
// const double Q11ee = m_collisions.Q11ee();
// const Eigen::ArrayXd& Q11ei = m_collisions.Q11ei();
// const Eigen::ArrayXd& Q12ei = m_collisions.Q12ei();
// const Eigen::ArrayXd& Q13ei = m_collisions.Q13ei();
// const Eigen::ArrayXd& Bstar = m_collisions.Bstij();
// const double Q22ee = m_collisions.Q22ee();
// double tauenergy=0.0;
//      double Vh0=0.0;
//      double collisionsum=0.0;
//    Vh0=std::sqrt(8*KB*Te/(PI*me));

//   collisionsum=(8.0/3.0)*(me/mh)*Vh0*nd*X[1]*Q11ei(1);

//   tauenergy=1/collisionsum;
//   return tauenergy;

// }

// //==============================================================================
// double Transport::tauEnergyBr()
// {
//  double nd  = m_thermo.numberDensity();
//   const double Te = m_thermo.Te();
//  const double me = m_thermo.speciesMw(0)/NA;
//   const double* const X = m_thermo.X();
//   const double B = m_thermo.getBField();
//   const double Th = m_thermo.T();

// double lam = 0.0;
//     double taue = 0.0;
//    double temporary = 0.0;
//     double TT=0.0;
//     double nne=0.0;
//     double mme=0.0;
//    double omegae= 0.0;

//     TT=Te*8.62E-5;
//    nne=nd*X[0]*1.0E-6;
//    temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//    if (TT>36.2){
//    temporary=25.25+std::log(TT)-0.5*std::log(nne);
//    }
//    taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//    taue=3.0/taue;
//    return taue;
//        }
// //==============================================================================
// //
// double Transport::etaohmbraginskii()
// {
//  double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//  const double me = m_thermo.speciesMw(0)/NA;
//   const double* const X = m_thermo.X();
//   const double B = m_thermo.getBField();
//   const double Th = m_thermo.T();
//          double lam = 0.0;
//         double taue = 0.0;
//          double temporary = 0.0;
//         double TT=0.0;
//        double nne=0.0;
//          double mme=0.0;
//                double omegae= 0.0;
//         double etaohm=0.0;
//        double factor= 0.0;
// factor=me/(QE*QE*nd*X[0]*MU0);

//        TT=Te*8.62E-5;
//        nne=nd*X[0]*1.0E-6;
//        temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//         if (TT>36.2){
//       temporary=25.25+std::log(TT)-0.5*std::log(nne);
//            }
//         taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//     etaohm=factor*(0.51)/(taue);


//      return etaohm;
// }

// //==============================================================================

// double Transport::transetaohmBr()
// {
//  double nd  = m_thermo.numberDensity();
//   const double Te = m_thermo.Te();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double* const X = m_thermo.X();
//    const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();

//    double transetaohmBr=0.0;
//    double lam = 0.0;
//    double taue = 0.0;
//    double temporary = 0.0;
//    double TT=0.0;
//    double nne=0.0;
//    double mme=0.0;
//    double omegae= 0.0;
//    double eta0=0.0;
//    double xx=0.0;
//    double delta=0.0;
//    omegae=QE*B/me;
//    TT=Te*8.62E-5;
//    nne=nd*X[0]*1.0E-6;
//    temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//    if (TT>36.2){
//    temporary=25.25+std::log(TT)-0.5*std::log(nne);
//    }
//    taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//    eta0=me/(QE*QE*nd*X[0]*taue*MU0);
//    xx=omegae*taue;
//    delta=xx*xx*xx*xx+14.8*xx*xx+3.77;

//    transetaohmBr=-eta0*xx*(1.7*xx*xx+0.78)/delta;
//    return transetaohmBr;
//   }
// //==============================================================================
// double Transport::perpetaohmBr()
// {
//  double nd  = m_thermo.numberDensity();
//   const double Te = m_thermo.Te();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double* const X = m_thermo.X();
//    const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();

//    double perpetaohmBr=0.0;
//    double lam = 0.0;
//    double taue = 0.0;
//    double temporary = 0.0;
//    double TT=0.0;
//    double nne=0.0;
//    double mme=0.0;
//    double omegae= 0.0;
//    double eta0=0.0;
//    double xx=0.0;
//    double delta=0.0;
//    omegae=QE*B/me;
//    TT=Te*8.62E-5;
//    nne=nd*X[0]*1.0E-6;
//    temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//    if (TT>36.2){
//    temporary=25.25+std::log(TT)-0.5*std::log(nne);
//    }
//    taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//    eta0=me/(QE*QE*nd*X[0]*taue*MU0);
//    xx=omegae*taue;
//    delta=xx*xx*xx*xx+14.8*xx*xx+3.77;

//    perpetaohmBr=eta0*(1-(6.42*xx*xx+1.84)/delta);
//    return perpetaohmBr;
//   }




// //==============================================================================
// double Transport::perpLambdaeBr()
// {
//  double nd  = m_thermo.numberDensity();
//  const double Te = m_thermo.Te();
//  const double me = m_thermo.speciesMw(0)/NA;
// const double* const X = m_thermo.X();
// const double B = m_thermo.getBField();
// const double Th = m_thermo.T();

// double perplambdaebr=0.0;
// double lam = 0.0;
// double taue = 0.0;
// double temporary = 0.0;
// double TT=0.0;
// double nne=0.0;
// double mme=0.0;
// double omegae= 0.0;
// double kappa0e=0.0;
// double xx=0.0;
// double delta=0.0;
// omegae=QE*B/me;
// TT=Te*8.62E-5;
// nne=nd*X[0]*1.0E-6;
// temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
// if (TT>36.2){
// temporary=25.25+std::log(TT)-0.5*std::log(nne);
// }
// taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
// kappa0e=nd*X[0]*KB*KB*Te*taue/me;
// xx=omegae*taue;
// delta=xx*xx*xx*xx+14.8*xx*xx+3.77;

// perplambdaebr=kappa0e*(4.66*xx*xx+11.9)/delta;

// return perplambdaebr;
// }

// //==============================================================================
// double Transport::transLambdaeBr()
// {
//  double nd  = m_thermo.numberDensity();
//   const double Te = m_thermo.Te();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double* const X = m_thermo.X();
//    const double B = m_thermo.getBField();
//    const double Th = m_thermo.T();

//    double translambdaebr=0.0;
//    double lam = 0.0;
//    double taue = 0.0;
//    double temporary = 0.0;
//    double TT=0.0;
//    double nne=0.0;
//    double mme=0.0;
//    double omegae= 0.0;
//    double kappa0e=0.0;
//    double xx=0.0;
//    double delta=0.0;
//    omegae=QE*B/me;
//    TT=Te*8.62E-5;
//    nne=nd*X[0]*1.0E-6;
//    temporary=23.4-0.5*std::log(nne)+1.5*std::log(TT);
//    if (TT>36.2){
//    temporary=25.25+std::log(TT)-0.5*std::log(nne);
//    }
//    taue=3.44E5*TT*std::sqrt(TT)/(temporary*nne);
//    kappa0e=nd*X[0]*KB*KB*Te*taue/me;
//    xx=omegae*taue;
//    delta=xx*xx*xx*xx+14.8*xx*xx+3.77;

//    translambdaebr=kappa0e*xx*(2.5*xx*xx+21.7)/delta;
//    return translambdaebr;
//   }



// //==============================================================================

// std::vector<double> Transport::parallelThermalDiffusionRatio2()
// {

// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
//      const Eigen::ArrayXd& lam01ei = m_collisions.L01ei();
//     const Eigen::ArrayXd& lam02ei = m_collisions.L02ei();
//       std::vector<double> kTi(ns);
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=fac*L(0,0);
// lam01=fac*L(0,1);
// lam02=fac*L(0,2);
// lam11=fac*L(1,1);
// lam12=fac*L(1,2);
// lam22=fac*L(2,2);
// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;
// kTi[0] = 2.5*Te/Th*X[0]*(lam01*lam22 - lam02*lam12)/(lam11*lam22 - lam12*lam12);
//     for (int i = 1; i < ns; ++i){
//         kTi[i] = -2.5*Te/Th*X[0]*(lam01ei[i]*lam22 - lam02ei[i]*lam12)/(lam11*lam22 - lam12*lam12);
//     }
//     return kTi;
// }

// //==============================================================================
// std::vector<double> Transport::perpThermalDiffusionRatio2()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
//      const Eigen::ArrayXd& lam01ei = m_collisions.L01ei();
//         const Eigen::ArrayXd& lam02ei = m_collisions.L02ei();
//       std::vector<double> kTi(ns);
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=fac*L(0,0);
// lam01=fac*L(0,1);
// lam02=fac*L(0,2);
// lam11=fac*L(1,1);
// lam12=fac*L(1,2);
// lam22=fac*L(2,2);
// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;

// kTi[0] = 2.5*Te/Th*X[0]*(lam01*(lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11)) - lam02*lam12*(lam11*lam22-lamB11*lamB22-lam12*lam12)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));

//     for (int i = 1; i < ns; ++i){
//         kTi[i] = -2.5*Te/Th*X[0]*(lam01ei[i]*(lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11)) - lam02ei[i]*lam12*(lam11*lam22-lamB11*lamB22-lam12*lam12)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
//     }
//     return kTi;
// }
// //==============================================================================
// std::vector<double> Transport::transverseThermalDiffusionRatio2()
// {
// const int ns = m_thermo.nGas();
//     const double Th = m_thermo.T();
//     const double Te = m_thermo.Te();
//     const double nd = m_thermo.numberDensity();
//     const double* const X = m_thermo.X();
//     const double me = m_thermo.speciesMw(0)/NA;
//     const double B = m_thermo.getBField();
//     const double P = m_thermo.P();
//      const Eigen::ArrayXd& lam01ei = m_collisions.L01ei();
//         const Eigen::ArrayXd& lam02ei = m_collisions.L02ei();
//       std::vector<double> kTi(ns);
// Matrix3d L = m_collisions.Lee<3>();
// const double fac = 75.*KB/(64.*X[0])*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));
// double fac2=(4.*X[0]/(25.*nd*KB));
// double fac3=25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//     double lam00 = 0.0;
//     double lam01 = 0.0;
//     double lam02 = 0.0;
//     double lam11 = 0.0;
//     double lam12 = 0.0;
//     double lam22 = 0.0;
//     double lamB00=0.0;
//     double lamB11=0.0;
//     double lamB22=0.0;

// lam00=fac*L(0,0);
// lam01=fac*L(0,1);
// lam02=fac*L(0,2);
// lam11=fac*L(1,1);
// lam12=fac*L(1,2);
// lam22=fac*L(2,2);
// lamB00 = QE*B/(KB*Te*fac3);
// lamB11 = 2.5*lamB00;
// lamB22 = 1.75*lamB11;
//     kTi[0] = -2.5*Te/Th*X[0]*(lam01*(lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11)) + lam02*lam12*(lam11*lamB22+lam22*lamB11)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));

//     for (int i = 1; i < ns; ++i){
//         kTi[i] = 2.5*Te/Th*X[0]*(lam01ei[i]*(lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11)) + lam02ei[i]*lam12*(lam11*lamB22+lam22*lamB11)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
//     }
//     return kTi;
// }

    } // namespace Transport
} // namespace Mutation

