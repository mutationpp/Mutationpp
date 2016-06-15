/**
 * @file Transport.cpp
 *
 * @brief Implements Transport class.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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
      mp_viscosity(NULL),
      mp_thermal_conductivity(NULL),
      mp_diffusion_matrix(NULL),
      mp_wrk1(NULL),
      mp_tag(NULL)
{ 
    // Load the viscosity calculator
    setViscosityAlgo(viscosity);
    
    // Load the thermal conductivity calculator
    setThermalConductivityAlgo(lambda);
    
    // Load the diffusion matrix calculator
    mp_diffusion_matrix = new Ramshaw(m_collisions);
    
    // Allocate work array storage
    mp_wrk1 = new double [m_thermo.nSpecies()*3];
    mp_wrk2 = mp_wrk1 + m_thermo.nSpecies();
    mp_wrk3 = mp_wrk2 + m_thermo.nSpecies();
    if(m_thermo.nEnergyEqns() > 1) {
        mp_tag  = new int [m_thermo.nEnergyEqns()*5];
        m_thermo.getTagModes(mp_tag);
    }
    
    //thermo.stateModel()->notifyOnUpdate(this);
}
    
//==============================================================================
    
Transport::~Transport()
{
    delete mp_viscosity;
    delete mp_thermal_conductivity;
    delete mp_diffusion_matrix;
    
    delete [] mp_wrk1;
    delete [] mp_tag;
}

//==============================================================================

void Transport::setViscosityAlgo(const std::string& algo)
{
    if (mp_viscosity != NULL) delete mp_viscosity;
    mp_viscosity = Factory<ViscosityAlgorithm>::create(algo, m_collisions);
}

//==============================================================================

double Transport::viscosity() { return mp_viscosity->viscosity(); }

//==============================================================================

void Transport::setThermalConductivityAlgo(const std::string& algo)
{
    if (mp_thermal_conductivity != NULL) delete mp_thermal_conductivity;
    mp_thermal_conductivity = Factory<ThermalConductivityAlgorithm>::create(
        algo, m_collisions);
}

//==============================================================================

double Transport::heavyThermalConductivity()
{
    return mp_thermal_conductivity->thermalConductivity();
}

//==============================================================================

void Transport::thermalDiffusionRatios(double* const p_k)
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

double Transport::electronThermalConductivity(int order)
{
    assert(order == 2 || order == 3);

	if (!m_thermo.hasElectrons())
        return 0.0;

	const double xe = m_thermo.X()[0];
    const double Te = m_thermo.Te();
    const double fac = 75.*KB/64.*std::sqrt(TWOPI*KB*Te/m_collisions.mass()(0));

    // 2nd order solution
    if (order == 2)
        return fac * xe / m_collisions.Lee<2>()(1,1);
    
    // 3rd order solution
    Matrix3d L = m_collisions.Lee<3>();
    return fac * xe * L(2,2) / (L(1,1)*L(2,2) - L(1,2)*L(1,2));
}

//==============================================================================

double Transport::internalThermalConductivity(double T)
{
    m_thermo.thermoDB()->cpint(T, mp_wrk1);
    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nSpecies()));
}


//==============================================================================

double Transport::rotationalThermalConductivity()
{
	m_thermo.speciesCpOverR(
        m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(), m_thermo.Tel(),
        NULL, NULL, mp_wrk1, NULL, NULL);

    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nSpecies()));
}

//==============================================================================

double Transport::vibrationalThermalConductivity()
{
	m_thermo.speciesCpOverR(
        m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(), m_thermo.Tel(),
        NULL, NULL, NULL, mp_wrk1, NULL);

    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nSpecies()));
}

//==============================================================================

double Transport::electronicThermalConductivity()
{
	m_thermo.speciesCpOverR(
        m_thermo.T(), m_thermo.Te(), m_thermo.Tr(), m_thermo.Tv(), m_thermo.Tel(),
        NULL, NULL, NULL, NULL, mp_wrk1);

    return euken(Map<const ArrayXd>(mp_wrk1, m_thermo.nSpecies()));
}

//==============================================================================

double Transport::reactiveThermalConductivity()
{
	// Compute dX_i/dT
    m_thermo.dXidT(mp_wrk1);
    
    // Compute the thermal diffusion ratios
    thermalDiffusionRatios(mp_wrk2);
    
    // Combine to get the driving forces
    for (int i = 0; i < m_thermo.nSpecies(); i++)
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
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        lambda -= mp_wrk1[i] / m_thermo.speciesMw(i) * mp_wrk2[i] * Y[i] * rho;
    
    return (RU * m_thermo.T() * lambda);
}

//==============================================================================

double Transport::soretThermalConductivity()
{
	// Compute dX_i/dT
    m_thermo.dXidT(mp_wrk1);
    
    // Compute the thermal diffusion ratios
    thermalDiffusionRatios(mp_wrk2);
    
    // Combine to get the driving forces
    for (int i = 0; i < m_thermo.nSpecies(); i++)
        mp_wrk1[i] += mp_wrk2[i] / m_thermo.T();
    
    // Compute the diffusion velocities
    double E;
    stefanMaxwell(mp_wrk1, mp_wrk1, E);
    
    double lambda = 0.0;
    for (int i = 0; i < m_thermo.nSpecies(); i++)
        lambda -= mp_wrk2[i]*mp_wrk1[i];
    
    return (m_thermo.P()*lambda);
}

//==============================================================================

void Transport::equilDiffFluxFacs(double* const p_F)
{
	// Get some state data
	const int ns = m_thermo.nSpecies();
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
	const int ns = m_thermo.nSpecies();
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
	const int ns = m_thermo.nSpecies();
	const double T  = m_thermo.T();

    // Compute the dXj/dT term
    m_thermo.dXidT(mp_wrk1);

	// Add thermal diffusion ratio term
	thermalDiffusionRatios(mp_wrk2);
	for (int i = 0; i < ns; ++i)
		mp_wrk1[i] += mp_wrk2[i]/T;

    // Compute the element averaged diffusion coefficients
	equilDiffFluxFacs(p_F);
}

//==============================================================================

void Transport::equilDiffFluxFacsZ(double* const p_F)
{
	const int ns = m_thermo.nSpecies();
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
    const double* const p_dp, double* const p_V, double& E)
{
	const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    
    // Electron offset
    const int k = m_thermo.hasElectrons() ? 1 : 0;
    const int nh = ns-k;

    Map<const ArrayXd> dp(p_dp, ns);
    Map<ArrayXd> V(p_V, ns), qi(mp_wrk3, ns);

    // Need to place a tolerance on X and Y
    const double tol = 1.0e-16;
    static ArrayXd X, Y; Y.resize(ns); // only gets resized the first time
    X = Map<const ArrayXd>(m_thermo.X(), ns).max(tol); // Place a tolerance on X
    m_thermo.convert<X_TO_Y>(&X[0], &Y[0]);

    // Get reference to binary diffusion coefficients
    const ArrayXd& nDij = m_collisions.nDij();

    // Compute mixture charge
    for (int i = 0; i < ns; ++i)
        qi[i] = m_thermo.speciesCharge(i);
    const double q = (X*qi).sum();

    // Compute kappas (store in V)
    V = (X*qi - q*Y)/(KB*Th);

    // Compute ambipolar electric field
    E = (k > 0 ? dp[0]/V[0] : 0.0);

    // Compute right hand side
    static VectorXd b; b.resize(nh);
    b.array() = -dp.tail(nh) + E*V.tail(nh);

    // Compute singular system matrix (lower triangular part)
    double fac;
    static MatrixXd G; G.resize(nh,nh);
    G.triangularView<Lower>() = MatrixXd::Zero(nh, nh);

    for (int j = 0, si = 1; j < nh; ++j, ++si) {
        for (int i = j+1; i < nh; ++i, ++si) {
            fac = X[i+k]*X[j+k]/nDij(si)*nd;
            G(i,i) += fac;
            G(j,j) += fac;
            G(i,j) = -fac;
        }
    }

    // Compute the constraints that are applied to the matrix
    fac = X[0]*qi[0];
    V.tail(nh) = Y.tail(nh);
    if (k > 0) V.tail(nh) -= X.tail(nh)*qi.tail(nh)*Y[0]/fac;

    // Add mass balance relation to make matrix nonsigular
    G.selfadjointView<Lower>().rankUpdate(
        V.matrix().tail(nh), nd/nDij.mean());

    static Eigen::LDLT<MatrixXd, Lower> ldlt;
    ldlt.compute(G);
    V.tail(nh) = ldlt.solve(b).array();

    // Compute electron diffusion velocity if need be
    if (k == 0)
        return;

    V[0] = -(X.tail(nh)*qi.tail(nh)*V.tail(nh)).sum() / fac;
}

//==============================================================================

double Transport::sigma(int order)
{
    assert(order == 1 || order == 2);

	if (!m_thermo.hasElectrons())
        return 0.0;

    const double xe = m_thermo.X()[0];
    const double me = m_collisions.mass()(0);
    const double Te = m_thermo.Te();
    const double fac = 3.*xe*QE*QE/(16.*KB*Te)*std::sqrt(TWOPI*KB*Te/me);
    
    // First order
    if (order == 1)
        return fac / m_collisions.Lee<1>()(0,0);
    
    // Second order
    Matrix2d L = m_collisions.Lee<2>();
    return fac / (L(0,0)-L(0,1)*L(0,1)/L(1,1));
}

template <>
Eigen::Vector3d Transport::electronThermalConductivityB<1>() {
    return Eigen::Vector3d::Zero();
}

//==============================================================================

//double Transport::meanFreePath()
//{
//	// Thermo properties
//    const int ns = m_thermo.nSpecies();
//    const int nh = m_thermo.nHeavy();
//    const int k  = ns-nh;
//    const double* const X = m_thermo.X();
//
//    // Get Q11 Collision integral
//    const ArrayXd& Q11 = m_collisions.Q11ij();
//
//    // Loop to sum over all heavy species
//    double sum = 0.0;
//    for (int i = 0; i < nh; ++i)
//        for (int j = i; j < nh; ++j)
//            sum += X[i]*X[j]*Q11(i,j);
//
//    return 1.0/(nd*sum);
//}
//
//
////==============================================================================
//
//double Transport::electronMeanFreePath()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons())
//        return 0.0;
//
//    // Thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//
//    // Get Q11 Collision integral
//    const MatrixXd& Q11 = mp_collisions->Q11(Th, Te, nd, X);
//
//    // Loop to sum over electrons and all species
//    double sum = 0.0;
//
//    for (int i = 0; i < ns; ++i)
//        sum += X[i]*X[0]*Q11(i);
//
//    return 1.0/(nd*sum);
//}
//
////==============================================================================
//
//double Transport::averageHeavyThermalSpeed()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	const int ns = m_thermo.nSpecies();
//    const int nh = m_thermo.nHeavy();
//    const double Th = m_thermo.T();
//    const double* const X = m_thermo.X();
//
//    // Loop to get average mass of heavy species
//    double ave_mw = 0.0;
//    for (int i = ns-nh; i < ns; ++i)
//        ave_mw += X[i] * m_thermo.speciesMw(i);
//
//    return sqrt(8.0*RU*Th/(PI*ave_mw));
//}
//
////==============================================================================
//
//double Transport::electronThermalSpeed()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	// Thermo properties
//    const double Te = m_thermo.Te();
//    const double Me = m_thermo.speciesMw(0);
//    return sqrt(8*RU*Te/(PI*Me));
//}
//
//
////==============================================================================
//
//double Transport::averageHeavyCollisionFreq()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	return meanFreePath()/averageHeavyThermalSpeed();
//}
//
////==============================================================================
//
//double Transport::electronHeavyCollisionFreq()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons())
//        return 0.0;
//    return electronThermalSpeed()/electronMeanFreePath();
//}
//
////==============================================================================
//
//double Transport::coulombMeanCollisionTime()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons())
//        return 0.0;
//
//    // Thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//
//    // Get Q11 Collision integral
//    const MatrixXd& Q11 = mp_collisions->Q11(Th, Te, nd, X);
//
//    // Loop to sum over electrons and all species
//    double sum = 0.0;
//
//    for (int i = 0; i < ns; ++i) {
//        sum = X[i]*Q11(i);
//    }
//
//    return (3.0/16.0)*1.0/(nd*sum);
//}
//
////==============================================================================
//
//double Transport::hallParameter()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons())
//        return 0.0;
//
//    // Thermo
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double B = m_thermo.getBField();
//
//    return QE*B*coulombMeanCollisionTime()/me;
//}
//
////==============================================================================
//
//double Transport::parallelDiffusionCoefficient()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	// Get thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double P = m_thermo.P();
//
//    // Get collision information
//    const MatrixXd& Q11   = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22   = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Cei   = mp_collisions->Cstei(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute the lambdas
//    double fac;
//    double fac1  = 0.0;
//    double fac2  = 0.0;
//    double lam00 = 0.0;
//    double lam01 = 0.0;
//    double lam02 = 0.0;
//    double lam11 = 0.0;
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//
//    double lamB00;
//    double lamB11;
//    double lamB22;
//
//    for (int i = 1; i < ns; ++i) {
//        fac1 = X[i]*Q11(i);
//        lam00 += fac1;
//        lam01 += fac1*(2.5-3.0*Cei(i));
//        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
//        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//
//    lam00 = fac*lam00;
//    lam01 = fac*lam01;
//    lam02 = fac*lam02;
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//
//    // First order
//    //return 1.0/lam00;
//
//    //Second order
//    //return 1/(lam00-lam01*lam01/lam11);
//
//    // Third order
//    return -(lam12*lam12-lam11*lam22)/(lam00*lam11*lam22-lam00*lam12*lam12-lam22*lam01*lam01+2.0*lam01*lam02*lam12-lam11*lam02*lam02);
//}
//
////==============================================================================
//
//double Transport::perpDiffusionCoefficient()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	// Get thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double P = m_thermo.P();
//    const double B = m_thermo.getBField();
//
//    // Get collision information
//    const MatrixXd& Q11   = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22   = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Cei   = mp_collisions->Cstei(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute the lambdas
//    double fac;
//    double fac1  = 0.0;
//    double fac2 = 0.0;
//    double lam00 = 0.0;
//    double lam01 = 0.0;
//    double lam02 = 0.0;
//    double lam11 = 0.0;
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//
//    double lamB00;
//    double lamB11;
//    double lamB22;
//
//    for (int i = 1; i < ns; ++i) {
//        fac1 = X[i]*Q11(i);
//        lam00 += fac1;
//        lam01 += fac1*(2.5-3.0*Cei(i));
//        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
//        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//    lam00 = fac*lam00;
//    lam01 = fac*lam01;
//    lam02 = fac*lam02;
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//
//    lamB00 = QE*B/(KB*Te*fac2);
//    lamB11 = 2.5*lamB00;
//    lamB22 = 1.75*lamB11;
//
//    // First order
//    //return lam00/(lam00*lam00+lamB00*lamB00);
//
//    // Second order
//   /* double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
//    double numerator = lam11*(lam00*lam11-lamB00*lamB11-lam01*lam01) + lamB11*(lamB00*lam11+lamB11*lam00);
//    return numerator/denominator;
//    */
//    // Third Order
//
//    double denominator = (lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02);
//    double numerator = -((lam12*lam12-lam11*lam22+lamB11*lamB22)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) - (lam11*lamB22+lam22*lamB11)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02));
//    return numerator/denominator;
//
//}
//
////==============================================================================
//
//double Transport::transverseDiffusionCoefficient()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	// Get thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double P = m_thermo.P();
//    const double B = m_thermo.getBField();
//
//    // Get collision information
//    const MatrixXd& Q11   = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22   = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar     = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Cei = mp_collisions->Cstei(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute the lambdas
//    double fac;
//    double fac1  = 0.0;
//    double fac2 = 0.0;
//    double lam00 = 0.0;
//    double lam01 = 0.0;
//    double lam02 = 0.0;
//    double lam11 = 0.0;
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//
//    double lamB00;
//    double lamB11;
//    double lamB22;
//
//    for (int i = 1; i < ns; ++i) {
//        fac1 = X[i]*Q11(i);
//        lam00 += fac1;
//        lam01 += fac1*(2.5-3.0*Cei(i));
//        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
//        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//    lam00 = fac*lam00;
//    lam01 = fac*lam01;
//    lam02 = fac*lam02;
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//
//    lamB00 = QE*B/(KB*Te*fac2);
//    lamB11 = 2.5*lamB00;
//    lamB22 = 1.75*lamB11;
//    // First order
//    //return -lamB00/(lam00*lam00+lamB00*lamB00);
//
//    // Second order
//   /* double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
//    double numerator = lamB11*(lam00*lam11-lamB00*lamB11-lam01*lam01) - lam11*(lamB00*lam11+lamB11*lam00);
//    return numerator/denominator;
//    */
//    // Third order
//
//    double denominator = (lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02);
//    double numerator = ((lam11*lamB22+lam22*lamB11)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam12*lam12-lam11*lam22+lamB11*lamB22)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02));
//    return numerator/denominator;
//
//}
////==============================================================================
//
//double Transport::parallelThermalDiffusionCoefficient()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    // Get thermodynamic properties
//    const int ns     = m_thermo.nSpecies();
//    const double Th  = m_thermo.T();
//    const double Te  = m_thermo.Te();
//    const double nd  = m_thermo.numberDensity();
//    const double me  = m_thermo.speciesMw(0)/NA;
//    const double *const X = m_thermo.X();
//
//    // Get collision integral information
//    const MatrixXd& Q11   = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22   = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar     = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute the lambdas
//    double fac;
//    double fac2 = 0.0;
//    double lam11 = 0.0;
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//
//    for (int i = 1; i < ns; ++i) {
//        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//    // Second order
//    return (lam22/(lam11*lam22-lam12*lam12));
//}
//
////==============================================================================
//
//double Transport::perpThermalDiffusionCoefficient()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    // Get thermodynamic properties
//    const int ns     = m_thermo.nSpecies();
//    const double Th  = m_thermo.T();
//    const double Te  = m_thermo.Te();
//    const double nd  = m_thermo.numberDensity();
//    const double me  = m_thermo.speciesMw(0)/NA;
//    const double *const X = m_thermo.X();
//    const double P   = m_thermo.P();
//    const double B = m_thermo.getBField();
//
//    // Get collision integral information
//    const MatrixXd& Q11   = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22   = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute the lambdas
//    double fac;
//    double fac2 = 0.0;
//    double lam11 = 0.0;
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//
//    double lamB00;
//    double lamB11;
//    double lamB22;
//
//    for (int i = 1; i < ns; ++i) {
//        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//    lamB00 = QE*B/(KB*Te*fac2);
//    lamB11 = 2.5*lamB00;
//    lamB22 = 1.75*lamB11;
//
//    // Second Order
//    double numerator = lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11);
//    double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);
//
//    return numerator/denominator;
//}
//
////==============================================================================
//
//double Transport::transverseThermalDiffusionCoefficient()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    // Get thermodynamic properties
//    const int ns     = m_thermo.nSpecies();
//    const double Th  = m_thermo.T();
//    const double Te  = m_thermo.Te();
//    const double nd  = m_thermo.numberDensity();
//    const double me  = m_thermo.speciesMw(0)/NA;
//    const double *const X = m_thermo.X();
//    const double P   = m_thermo.P();
//    const double B = m_thermo.getBField();
//
//    // Get collision integral information
//    const MatrixXd& Q11   = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22   = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute the lambdas
//    double fac;
//    double fac2 = 0.0;
//    double lam11 = 0.0;
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//
//    double lamB00;
//    double lamB11;
//    double lamB22;
//
//    for (int i = 1; i < ns; ++i) {
//        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//    lamB00 = QE*B/(KB*Te*fac2);
//    lamB11 = 2.5*lamB00;
//    lamB22 = 1.75*lamB11;
//
//    // Second Order
//    double numerator = lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11);
//    double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);
//
//    return numerator/denominator;
//}
////==============================================================================
//
//double Transport::sigmaParallel()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double* const X = m_thermo.X();
//    const double nd = m_thermo.numberDensity();
//    const double P = m_thermo.P();
//
//
//    return 4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te)*(nd*KB*Te/P)* parallelDiffusionCoefficient();
//
//}
////==============================================================================
//
//double Transport::sigmaPerpendicular()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double* const X = m_thermo.X();
//    const double nd = m_thermo.numberDensity();
//    const double P = m_thermo.P();
//
//    return 4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te)*(nd*KB*Te/P)* perpDiffusionCoefficient();
//
//}
//
////==============================================================================
//
//double Transport::sigmaTransverse()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double* const X = m_thermo.X();
//    const double nd = m_thermo.numberDensity();
//    const double P = m_thermo.P();
//
//    return -4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te)*(nd*KB*Te/P)* transverseDiffusionCoefficient();
//
//}
//
////==============================================================================
//
//double Transport::parallelElectronThermalConductivity()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    const double* const X = m_thermo.X();
//    const double nd  = m_thermo.numberDensity();
//    const double Th  = m_thermo.T();
//    const double Te  = m_thermo.Te();
//
//    double Kt_parallel = parallelThermalDiffusionCoefficient();
//    return (X[0]*X[0])*Kt_parallel;
//}
//
////==============================================================================
//
//double Transport::perpElectronThermalConductivity()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    const double* const X = m_thermo.X();
//    const double nd  = m_thermo.numberDensity();
//    const double Th  = m_thermo.T();
//    const double Te  = m_thermo.Te();
//
//
//    return (X[0]*X[0])*perpThermalDiffusionCoefficient();
//}
//
////==============================================================================
//
//double Transport::transverseElectronThermalConductivity()
//{
//	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)
//
//    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
//        return 0.0;
//
//    const double* const X = m_thermo.X();
//    const double nd  = m_thermo.numberDensity();
//    const double Th  = m_thermo.T();
//    const double Te  = m_thermo.Te();
//
//    return -(X[0]*X[0])*transverseThermalDiffusionCoefficient();
//}
//
////==============================================================================
//
//double Transport::ratioSigmaPerpPar()
//{
//    return sigmaPerpendicular()/sigmaParallel();
//}
////==============================================================================
//
//double Transport::ratioSigmaTransPar()
//{
//    return sigmaTransverse()/sigmaParallel();
//}
////==============================================================================
//
//double Transport::ratioLambdaPerpPar()
//{
//    return perpElectronThermalConductivity()/parallelElectronThermalConductivity();
//}
////==============================================================================
//
//double Transport::ratioLambdaTransPar()
//{
//    return transverseElectronThermalConductivity()/parallelElectronThermalConductivity();
//}
////==============================================================================
//
//std::vector<double> Transport::parallelThermalDiffusionRatio()
//{
//    double* species_values = new double [m_thermo.nSpecies()];
//
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//    const double me = m_thermo.speciesMw(0)/NA;
//
//    const MatrixXd& Q11 = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22 = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar   = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Cei = mp_collisions->Cstei(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute lambdas
//    double lam00 = 0.0;
//    double lam01 = 0.0;
//    double lam11 = 0.0;
//    double lam02 = 0.0;
//    std::vector<double> lam01ei(ns);
//    std::vector<double> lam02ei(ns);
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//    double fac   = 0.0;
//
//
//    for (int i = 1; i < ns; ++i) {
//        fac = X[i]*Q11(i);
//        lam00 += fac;
//        lam01 += fac*(2.5-3.0*Cei(i));
//        lam11 += fac*(25.0/4.0-3.0*Bstar(i));
//        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//    for (int i = 0; i < ns; ++i) {
//    lam01ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(2.5*Q11(i)-3.0*Q12ei(i));
//    lam02ei[i] =64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(35.0/8.0*Q11(i) - 21.0/2.0*Q12ei(i) + 6.0*Q13ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    lam00 = fac*lam00;
//    lam01 = fac*lam01;
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam02 = fac*lam02;
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//    std::vector<double> kTi(ns);
//    // Second Order
//   /* kTi[0] = 2.5*Te/Th*X[0]*lam01/lam11;
//    for (int i = 1; i < ns; ++i){
//        kTi[i] = -2.5*Te/Th*X[0]*lam01ei[i]/lam11;
//    }
//    */
//
//    // Third Order
//    kTi[0] = 2.5*Te/Th*X[0]*(lam01*lam22 - lam02*lam12)/(lam11*lam22 - lam12*lam12);
//    for (int i = 1; i < ns; ++i){
//        kTi[i] = -2.5*Te/Th*X[0]*(lam01ei[i]*lam22 - lam02ei[i]*lam12)/(lam11*lam22 - lam12*lam12);
//    }
//    return kTi;
//
//}
////==============================================================================
//std::vector<double> Transport::perpThermalDiffusionRatio()
//{
//    // Get thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double B = m_thermo.getBField();
//
//    // Get collision integrals
//    const MatrixXd& Q11 = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22 = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar   = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Cei = mp_collisions->Cstei(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute lambdas
//    double lam00 = 0.0;
//    double lam01 = 0.0;
//    double lam11 = 0.0;
//    double lam02 = 0.0;
//    std::vector<double> lam01ei(ns);
//    std::vector<double> lam02ei(ns);
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//    double fac   = 0.0;
//    double fac2 = 0.0;
//
//    double lamB00 = 0.0;
//    double lamB11 = 0.0;
//    double lamB22 = 0.0;
//
//    for (int i = 1; i < ns; ++i) {
//        fac = X[i]*Q11(i);
//        lam00 += fac;
//        lam01 += fac*(2.5-3.0*Cei(i));
//        lam11 += fac*(25.0/4.0-3.0*Bstar(i));
//        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//
//    for (int i = 0; i < ns; ++i) {
//        lam01ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(2.5*Q11(i)-3.0*Q12ei(i));
//        lam02ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(35.0/8.0*Q11(i) - 21.0/2.0*Q12ei(i) + 6.0*Q13ei(i));
//    }
//
//    lam00 = fac*lam00;
//    lam01 = fac*lam01;
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam02 = fac*lam02;
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//    lamB00 = QE*B/(KB*Te*fac2);
//    lamB11 = 2.5*lamB00;
//    lamB22 = 1.75*lamB11;
//
//    std::vector<double> kTi(ns);
//
//    // Third Order
//
//    kTi[0] = 2.5*Te/Th*X[0]*(lam01*(lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11)) - lam02*lam12*(lam11*lam22-lamB11*lamB22-lam12*lam12)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
//
//    for (int i = 1; i < ns; ++i){
//        kTi[i] = -2.5*Te/Th*X[0]*(lam01ei[i]*(lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11)) - lam02ei[i]*lam12*(lam11*lam22-lamB11*lamB22-lam12*lam12)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
//    }
//    return kTi;
//
//}
////==============================================================================
//
//std::vector<double> Transport::transverseThermalDiffusionRatio()
//{
//    // Get thermo properties
//    const int ns = m_thermo.nSpecies();
//    const double Th = m_thermo.T();
//    const double Te = m_thermo.Te();
//    const double nd = m_thermo.numberDensity();
//    const double* const X = m_thermo.X();
//    const double me = m_thermo.speciesMw(0)/NA;
//    const double B = m_thermo.getBField();
//
//    // Get collision integrals
//    const MatrixXd& Q11 = mp_collisions->Q11(Th, Te, nd, X);
//    const MatrixXd& Q22 = mp_collisions->Q22(Th, Te, nd, X);
//    const MatrixXd& Bstar   = mp_collisions->Bstar(Th, Te, nd, X);
//    const VectorXd& Cei = mp_collisions->Cstei(Th, Te, nd, X);
//    const VectorXd& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
//    const VectorXd& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
//    const VectorXd& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
//    const VectorXd& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
//    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
//    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
//
//    // Compute lambdas
//    double lam00 = 0.0;
//    double lam01 = 0.0;
//    double lam11 = 0.0;
//    double lam02 = 0.0;
//    std::vector<double> lam01ei(ns);
//    std::vector<double> lam02ei(ns);
//    double lam12 = 0.0;
//    double lam22 = 0.0;
//    double fac   = 0.0;
//    double fac2 = 0.0;
//
//    double lamB00 = 0.0;
//    double lamB11 = 0.0;
//    double lamB22 = 0.0;
//
//    for (int i = 1; i < ns; ++i) {
//        fac = X[i]*Q11(i);
//        lam00 += fac;
//        lam01 += fac*(2.5-3.0*Cei(i));
//        lam11 += fac*(25.0/4.0-3.0*Bstar(i));
//        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
//        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
//                       Q14ei(i));
//        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
//                       210.0*Q14ei(i)+90.0*Q15ei(i));
//    }
//
//    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
//    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
//
//    for (int i = 0; i < ns; ++i) {
//        lam01ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(2.5*Q11(i)-3.0*Q12ei(i));
//        lam02ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(35.0/8.0*Q11(i) - 21.0/2.0*Q12ei(i) + 6.0*Q13ei(i));
//    }
//
//    lam00 = fac*lam00;
//    lam01 = fac*lam01;
//    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
//    lam02 = fac*lam02;
//    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
//    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
//
//    lamB00 = QE*B/(KB*Te*fac2);
//    lamB11 = 2.5*lamB00;
//    lamB22 = 1.75*lamB11;
//
//    std::vector<double> kTi(ns);
//
//    // Third Order
//
//    kTi[0] = -2.5*Te/Th*X[0]*(lam01*(lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11)) + lam02*lam12*(lam11*lamB22+lam22*lamB11)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
//
//    for (int i = 1; i < ns; ++i){
//        kTi[i] = 2.5*Te/Th*X[0]*(lam01ei[i]*(lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11)) + lam02ei[i]*lam12*(lam11*lamB22+lam22*lamB11)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
//    }
//    return kTi;
//
//}
////==============================================================================
//std::vector<double> Transport::ratiokTPerpPar()
//{
//    const int ns = m_thermo.nSpecies();
//
//    std::vector<double> ratio(ns);
//    for (int i = 0; i < ns; ++i){
//        ratio[i] = perpThermalDiffusionRatio()[i]/parallelThermalDiffusionRatio()[i];
//    }
//    return ratio;
//}
////==============================================================================
//std::vector<double> Transport::ratiokTTransPar()
//{
//    const int ns = m_thermo.nSpecies();
//
//    std::vector<double> ratio(ns);
//    for (int i = 0; i < ns; ++i){
//        ratio[i] = transverseThermalDiffusionRatio()[i]/parallelThermalDiffusionRatio()[i];
//    }
//    return ratio;
//}
////==============================================================================

    } // namespace Transport
} // namespace Mutation

