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

#include "Transport.h"
#include "Constants.h"
#include "Numerics.h"

#include <iostream>

using namespace Mutation::Numerics;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace Transport {

using Mutation::Thermodynamics::Thermodynamics;

//==============================================================================

#define ERROR_IF_INTEGRALS_ARE_NOT_LOADED(__RET__)\
if (mp_collisions == NULL) {\
	cout << "Error! Trying to use transport without loading collision integrals!!" << endl;\
	return __RET__;\
}

//==============================================================================

Transport::Transport(
    Thermodynamics& thermo, const std::string& viscosity, const std::string& lambda, const bool load_data)
    : m_thermo(thermo),
      mp_collisions(NULL),
      mp_viscosity(NULL),
      mp_thermal_conductivity(NULL),
      mp_diffusion_matrix(NULL),
      mp_wrk1(NULL),
      mp_tag(NULL)
{ 
    if (!load_data)
        return;

    // Load the collision integral data
    mp_collisions = new CollisionDB(thermo);

    // Load the viscosity calculator
    mp_viscosity = 
        Config::Factory<ViscosityAlgorithm>::create(viscosity, *mp_collisions);
    
    // Load the thermal conductivity calculator
    mp_thermal_conductivity = 
        Config::Factory<ThermalConductivityAlgorithm>::create(
            lambda, *mp_collisions);
    
    // Load the diffusion matrix calculator
    mp_diffusion_matrix =
        new Ramshaw(thermo, *mp_collisions);
    
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
    delete mp_collisions;
    delete mp_viscosity;
    delete mp_thermal_conductivity;
    delete mp_diffusion_matrix;
    
    delete [] mp_wrk1;
    delete [] mp_tag;
}

//==============================================================================

void Transport::omega11ii(double* const p_omega)
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

    const double ns  = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double *const X = m_thermo.X();
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    
    for (int i = 0; i < ns; ++i)
        p_omega[i] = Q11(i,i);
}

//==============================================================================

void Transport::omega22ii(double* const p_omega)
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

    const double ns  = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double *const X = m_thermo.X();
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    
    for (int i = 0; i < ns; ++i)
        p_omega[i] = Q22(i,i);
}

//==============================================================================

void Transport::frozenThermalConductivityVector(double* const p_lambda)
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

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

double Transport::electronThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    // Get thermodynamic properties
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double me  = m_thermo.speciesMw(0)/NA;
    const double *const X = m_thermo.X();
    
    // Get collision integral information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& B     = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);
    
    // Compute the lambdas
    double fac;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;
    
    for (int i = 1; i < ns; ++i) {
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*B(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
            Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
            210.0*Q14ei(i)+90.0*Q15ei(i));
    }
    
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));
    
    if (lam11 <= 0.0) {
        cout << "Negative electron thermal conductivity!" << endl;
        cout << "electron mole fraction: " << X[0] << endl;
        cout << "heavy temperature: " << Th << endl;
        cout << "electron temperature: " << Te << endl;
        cout << "Q22(0): " << Q22(0) << endl;
        cout << "species, 6.25 - 3.0*B(i) > 0, Q11(i) > 0" << endl;
        for (int i = 1; i < ns; ++i)
            cout << m_thermo.speciesName(i) << ", " << 6.25 - 3.0*B(i) << ", " << Q11(i) << endl;
    }
    assert(lam11 > 0.0);
    //assert(lam22 > 0.0);

    //if (lam11*lam22 <= lam12*lam12)
    //    cout << lam11*lam22 << " " << lam12*lam12 << " " << X[0] << endl;
    //assert(lam11*lam22 > lam12*lam12);

    // 2nd order solution
    //return (X[0]*X[0]/std::max(lam11, 1.0e-20));
    
    // 3rd order solution
    double denom = (lam11*lam22 > lam12*lam12 ? lam11*lam22-lam12*lam12 : 1.0e-30);
    return std::max(0.0, (X[0]*X[0]*lam22/denom));
}

//==============================================================================

double Transport::internalThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, mp_wrk1, mp_wrk2, mp_wrk3);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * (mp_wrk1[i] + mp_wrk2[i] + mp_wrk3[i]) / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::rotationalThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, mp_wrk1, NULL, NULL);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * mp_wrk1[i] / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::vibrationalThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, NULL, mp_wrk1, NULL);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * mp_wrk1[i] / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================

double Transport::electronicThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    const int ns     = m_thermo.nSpecies();
    const int nh     = m_thermo.nHeavy();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double Tr  = m_thermo.Tr();
    const double Tv  = m_thermo.Tv();
    const double Tel = m_thermo.Tel();
    const double nd  = m_thermo.numberDensity(); 
    const double *const X = m_thermo.X();
    
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, X);
    
    m_thermo.speciesCpOverR(
        Th, Te, Tr, Tv, Tel, NULL, NULL, NULL, NULL, mp_wrk1);
    
    double lambda = 0.0;
    for (int i = ns-nh; i < ns; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ns; ++j)
            sum += X[j] / nDij(i,j);
        
        lambda += X[i] * mp_wrk1[i] / sum;
    }
    
    return (lambda * KB);
}

//==============================================================================
/*double Transport::reactiveThermalConductivity()
{
    ///
    /// @todo This function fails if there is no electric field.  Does not take
    /// into account what happens when there are no electrons.
    ///

    // Compute dX_i/dT
    m_thermo.dXidT(m_thermo.T(), m_thermo.P(), m_thermo.X(), mp_wrk1);
    
    // Compute the multicomponent diffusion coefficient matrix
    const int ns = m_thermo.nSpecies();
    const RealMatrix& Dij = diffusionMatrix();
    const double* const X = m_thermo.X();
    const double* const Y = m_thermo.Y();
    
    // Store the species charges and determine the mixture charge
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }
    
    for (int i = 0; i < ns; ++i) {
        mp_wrk2[i] = 0.0;
        for (int j = 0; j < ns; ++j) {
            mp_wrk2[i] += Dij(i,j)*(mp_wrk1[j]-(X[j]*mp_wrk3[j]-Y[j]*q)/
                (X[0]*mp_wrk3[0]-Y[0]*q)*mp_wrk1[0]);
        }
    }
    
    // Compute the species enthalpies per unit mass
    m_thermo.speciesHOverRT(mp_wrk1);
    
    double lambda = 0.0;
    for (int i = 0; i < ns; ++i)
        lambda += mp_wrk1[i] * mp_wrk2[i] * X[i];
    
    return (m_thermo.P()*lambda);
}*/

double Transport::reactiveThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

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

double Transport::soretThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

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

void Transport::averageDiffusionCoeffs(double *const p_Di)
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const p_X = m_thermo.X();
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, p_X);
    
    for (int i = 0; i < ns; ++i)
        p_Di[i] = 0.0;
    
    for (int i = 0, index = 1; i < ns; ++i, ++index) {
        for (int j = i + 1; j < ns; ++j, ++index) {
            p_Di[i] += p_X[j] / nDij(index);
            p_Di[j] += p_X[i] / nDij(index);
        }
    }
    
    for (int i = 0; i < ns; ++i)
        p_Di[i] = (1.0 - p_X[i]) / (p_Di[i] * nd);
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

	const RealMatrix& Dij = diffusionMatrix();
	const RealMatrix& nu  = m_thermo.elementMatrix();

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
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

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
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

	const int ns = m_thermo.nSpecies();
	const double T  = m_thermo.T();
	const double nd = m_thermo.numberDensity();
	const double* const p_X = m_thermo.X();

    // Compute the dXj/dT term
    m_thermo.dXidT(mp_wrk1);

	// Add thermal diffusion ratio term
	mp_thermal_conductivity->thermalDiffusionRatios(T, T, nd, p_X, mp_wrk2);
	for (int i = 0; i < ns; ++i)
		mp_wrk1[i] += mp_wrk2[i]/T;

    // Compute the element averaged diffusion coefficients
	equilDiffFluxFacs(p_F);
}

//==============================================================================

void Transport::equilDiffFluxFacsZ(double* const p_F)
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

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
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED()

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    
    // Need to place a tolerance on X and Y
    const double tol = 1.0e-16;
    static std::vector<double> X(ns);
    static std::vector<double> Y(ns);
    for (int i = 0; i < ns; ++i)
        X[i] = std::max(tol, m_thermo.X()[i]);
    m_thermo.convert<X_TO_Y>(&X[0], &Y[0]);

    // Get reference to binary diffusion coefficients
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, &X[0]);

    // Electron offset
    const int k = m_thermo.hasElectrons() ? 1 : 0;

    // Compute mixture charge
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }

    // Compute kappas
    const double one_over_kbt = 1.0/(KB*Th);
    for (int i = 0; i < ns; ++i)
        p_V[i] = (X[i]*mp_wrk3[i] - Y[i]*q)*one_over_kbt;

    // Compute ambipolar electric field
    E = (k > 0 ? p_dp[0]/p_V[0] : 0.0);

    // Compute right hand side
    static RealVector b(ns-k);
    for (int i = k; i < ns; ++i)
        b(i-k) = -p_dp[i] + p_V[i]*E;

    // Compute singular system matrix
    double fac;
    static RealSymMat G(ns-k);
    for (int i = k; i < ns; ++i)
        G(i-k,i-k) = 0.0;
    for (int i = k; i < ns; ++i) {
        for (int j = i+1; j < ns; ++j) {
            fac = X[i]*X[j]/nDij(i,j)*nd;
            G(i-k,i-k) += fac;
            G(j-k,j-k) += fac;
            G(i-k,j-k) = -fac;
        }
    }

    // Compute average binary diffusion coefficient
    double a = 0.0;
    for (int i = 0; i < nDij.size(); ++i)
        a += nDij(i);
    a /= ((double)nDij.size()*nd);
    a = 1.0/a;

    // Compute the constraints that are applied to the matrix
    fac = (k > 0 ? Y[0]/(X[0]*mp_wrk3[0]) : 0.0);
    for (int i = k; i < ns; ++i)
        p_V[i] = Y[i] - X[i]*mp_wrk3[i]*fac;

    // Add mass balance relation to make make matrix nonsigular
    for (int i = k; i < ns; ++i)
        for (int j = i; j < ns; ++j)
            G(i-k,j-k) += a*p_V[i]*p_V[j];

    // Solve for the velocities
    static RealVector V(ns-k);
    static LDLT<double> ldlt;
    ldlt.setMatrix(G);
    ldlt.solve(V, b);

    // Copy to output vector
    for (int i = k; i < ns; ++i)
        p_V[i] = V(i-k);

    // Compute electron diffusion velocity if need be
    if (k == 0)
        return;

    p_V[0] = 0.0;
    for (int i = k; i < ns; ++i)
        p_V[0] += X[i]*mp_wrk3[i]*p_V[i];
    p_V[0] /= -mp_wrk3[0]*X[0];
}

// For now assume that we have ions and solve the full system in thermal
// nonequilibrium with ambipolar assumption
/*void Transport::stefanMaxwell(
    const double* const p_dp, double* const p_V, double& E)
{
    // Determine some constants
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();

    // Place a tolerance on X and Y
    const double tol = 1.0e-16;
    static std::vector<double> xy(2*ns);
    double *X = &xy[0], *Y = &xy[ns];
    double sum = 0.0;
    for (int i = 0; i < ns; ++i) {
        X[i] = tol + m_thermo.X()[i];
        sum += X[i];
    }
    for (int i = 0; i < ns; ++i)
        X[i] /= sum;
    m_thermo.convert<X_TO_Y>(X, Y);

    // Get reference to binary diffusion coefficients
    const RealSymMat& nDij = mp_collisions->nDij(Th, Te, nd, &X[0]);

    // Compute mixture charge (store species' charges in mp_wrk3 work array)
    double q = 0.0;
    for (int i = 0; i < ns; ++i) {
        mp_wrk3[i] = m_thermo.speciesCharge(i);
        q += mp_wrk3[i] * X[i];
    }

    // Compute kappas (stored in p_V to avoid creating another array)
    // also compute the 2-norm of the kappa vector
    const double one_over_kbt = 1.0/(KB*Th);
    double s = 0.0;
    for (int i = 0; i < ns; ++i) {
        p_V[i] = (X[i]*mp_wrk3[i] - Y[i]*q)*one_over_kbt;
        s += p_V[i]*p_V[i];
    }
    s = std::sqrt(s);

    // Compute the RHS vector
    static RealVector b(ns+1);
    for (int i = 0; i < ns; ++i)
        b(i) = -p_dp[i];
    b(0) *= Th/Te;
    b(ns) = 0.0;
    //cout << b << endl;
    // Compute system matrix
    static RealMatrix G(ns+1,ns+1);

    // - electron contribution
    double fac1 = Te/Th;
    double fac2 = fac1*X[0]*nd;
    G(0,0) = 0.0;
    for (int j = 1; j < ns; ++j) {
        G(0,j) =  -fac2*X[j]/nDij(0,j);
        G(0,0) += G(0,j);
        G(j,0) =  G(0,j);
        G(j,j) = -fac1*G(j,0);
    }
    G(0,0) /= -fac1;
    G(0,ns) = -p_V[0]/(fac1*s);

    // - heavy contribution
    for (int i = 1; i < ns; ++i) {
        for (int j = i+1; j < ns; ++j) {
            G(i,j) =  -X[i]*X[j]/nDij(i,j)*nd;
            G(i,i) -= G(i,j);
            G(j,i) =  G(i,j);
            G(j,j) -= G(j,i);
        }
        G(i,ns) = -p_V[i]/s;
    }

    // - ambipolar constraint
    for (int j = 0; j < ns; ++j)
        G(ns,j) = -p_V[j]/s;
    G(ns,ns) = 0.0;

    //cout << G << endl;

    // Finally solve the system for Vi and E
    static RealVector x(ns+1);
    std::pair<int, double> ret = Numerics::gmres(
        G, x, b, Numerics::DiagonalPreconditioner<double>(G));

    // Compute mass constraint projector
    static RealVector R(ns,1.0);
    R(0) /= fac1;
    double r = 0.0;
    for (int i = 0; i < ns; ++i)
        r += R(i)*Y[i];
    double p = 0.0;
    for (int i = 0; i < ns; ++i)
        p += x(i)*Y[i];
    p /= r;

    // Retrieve the solution
    for (int i = 0; i < ns; ++i)
        p_V[i] = x(i) - p*R(i);
    E = b(ns)/s;
}*/

//==============================================================================

double Transport::sigma() 
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30) 
        return 0.0;

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;

    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& B   = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei = mp_collisions->Cstei(Th, Te, nd, X);
    
    // Compute lambdas
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = 0.0;
    double fac   = 0.0;
    
    for (int i = 1; i < ns; ++i) {
        fac = X[i]*Q11(i);
        lam00 += fac;
        lam01 += fac*(2.5-3.0*Cei(i));
        lam11 += fac*(25.0/4.0-3.0*B(i));
    }
    
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    
    // First order
    //return (4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te*lam00));
    
    // Second order
    return (4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te*(lam00-lam01*lam01/lam11)));
}

//==============================================================================

double Transport::meanFreePath()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	// Thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;

    // Get Q11 Collision integral
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);

    // Loop to sum over electrons and all species
    double sum = 0.0;

    for (int i = 0; i < ns; ++i)
        for (int j = 0; j < ns; ++j)
            sum += X[i]*X[j]*Q11(i,j);

    return 1.0/(nd*sum);
}


//==============================================================================

double Transport::electronMeanFreePath()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons())
        return 0.0;

    // Thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();

    // Get Q11 Collision integral
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);

    // Loop to sum over electrons and all species
    double sum = 0.0;

    for (int i = 0; i < ns; ++i)
        sum += X[i]*X[0]*Q11(i);

    return 1.0/(nd*sum);
}

//==============================================================================

double Transport::averageHeavyThermalSpeed()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	const int ns = m_thermo.nSpecies();
    const int nh = m_thermo.nHeavy();
    const double Th = m_thermo.T();
    const double* const X = m_thermo.X();

    // Loop to get average mass of heavy species
    double ave_mw = 0.0;
    for (int i = ns-nh; i < ns; ++i)
        ave_mw += X[i] * m_thermo.speciesMw(i);

    return sqrt(8.0*RU*Th/(PI*ave_mw));
}

//==============================================================================

double Transport::electronThermalSpeed()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	// Thermo properties
    const double Te = m_thermo.Te();
    const double Me = m_thermo.speciesMw(0);
    return sqrt(8*RU*Te/(PI*Me));
}


//==============================================================================

double Transport::averageHeavyCollisionFreq()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	return meanFreePath()/averageHeavyThermalSpeed();
}

//==============================================================================

double Transport::electronHeavyCollisionFreq()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons())
        return 0.0;
    return electronThermalSpeed()/electronMeanFreePath();
}

//==============================================================================

double Transport::coulombMeanCollisionTime()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons())
        return 0.0;

    // Thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();

    // Get Q11 Collision integral
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);

    // Loop to sum over electrons and all species
    double sum = 0.0;

    for (int i = 0; i < ns; ++i) {
        sum = X[i]*Q11(i);
    }

    return (3.0/16.0)*1.0/(nd*sum);
}

//==============================================================================

double Transport::hallParameter()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons())
        return 0.0;

    // Thermo
    const double me = m_thermo.speciesMw(0)/NA;
    const double B = m_thermo.getBField();

    return QE*B*coulombMeanCollisionTime()/me;
}

//==============================================================================

double Transport::parallelDiffusionCoefficient()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	// Get thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;
    const double P = m_thermo.P();

    // Get collision information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei   = mp_collisions->Cstei(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute the lambdas
    double fac;
    double fac1  = 0.0;
    double fac2  = 0.0;
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam02 = 0.0;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;

    double lamB00;
    double lamB11;
    double lamB22;

    for (int i = 1; i < ns; ++i) {
        fac1 = X[i]*Q11(i);
        lam00 += fac1;
        lam01 += fac1*(2.5-3.0*Cei(i));
        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];

    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam02 = fac*lam02;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));


    // First order
    //return 1.0/lam00;

    //Second order
    //return 1/(lam00-lam01*lam01/lam11);

    // Third order
    return -(lam12*lam12-lam11*lam22)/(lam00*lam11*lam22-lam00*lam12*lam12-lam22*lam01*lam01+2.0*lam01*lam02*lam12-lam11*lam02*lam02);
}

//==============================================================================

double Transport::perpDiffusionCoefficient()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	// Get thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;
    const double P = m_thermo.P();
    const double B = m_thermo.getBField();

    // Get collision information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei   = mp_collisions->Cstei(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute the lambdas
    double fac;
    double fac1  = 0.0;
    double fac2 = 0.0;
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam02 = 0.0;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;

    double lamB00;
    double lamB11;
    double lamB22;

    for (int i = 1; i < ns; ++i) {
        fac1 = X[i]*Q11(i);
        lam00 += fac1;
        lam01 += fac1*(2.5-3.0*Cei(i));
        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam02 = fac*lam02;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));


    lamB00 = QE*B/(KB*Te*fac2);
    lamB11 = 2.5*lamB00;
    lamB22 = 1.75*lamB11;

    // First order
    //return lam00/(lam00*lam00+lamB00*lamB00);

    // Second order
   /* double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
    double numerator = lam11*(lam00*lam11-lamB00*lamB11-lam01*lam01) + lamB11*(lamB00*lam11+lamB11*lam00);
    return numerator/denominator;
    */
    // Third Order

    double denominator = (lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02);
    double numerator = -((lam12*lam12-lam11*lam22+lamB11*lamB22)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) - (lam11*lamB22+lam22*lamB11)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02));
    return numerator/denominator;

}

//==============================================================================

double Transport::transverseDiffusionCoefficient()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	// Get thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;
    const double P = m_thermo.P();
    const double B = m_thermo.getBField();

    // Get collision information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar     = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei = mp_collisions->Cstei(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute the lambdas
    double fac;
    double fac1  = 0.0;
    double fac2 = 0.0;
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam02 = 0.0;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;

    double lamB00;
    double lamB11;
    double lamB22;

    for (int i = 1; i < ns; ++i) {
        fac1 = X[i]*Q11(i);
        lam00 += fac1;
        lam01 += fac1*(2.5-3.0*Cei(i));
        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }
    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam02 = fac*lam02;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));


    lamB00 = QE*B/(KB*Te*fac2);
    lamB11 = 2.5*lamB00;
    lamB22 = 1.75*lamB11;
    // First order
    //return -lamB00/(lam00*lam00+lamB00*lamB00);

    // Second order
   /* double denominator = (lam00*lam11-lamB00*lamB11-lam01*lam01)*(lam00*lam11-lamB00*lamB11-lam01*lam01) + (lamB00*lam11+lamB11*lam00)*(lamB00*lam11+lamB11*lam00);
    double numerator = lamB11*(lam00*lam11-lamB00*lamB11-lam01*lam01) - lam11*(lamB00*lam11+lamB11*lam00);
    return numerator/denominator;
    */
    // Third order

    double denominator = (lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02);
    double numerator = ((lam11*lamB22+lam22*lamB11)*(lam00*lam11*lam22-lam22*lamB00*lamB11-lamB22*(lam00*lamB11+lam11*lamB00)-lam00*lam12*lam12-lam01*lam01*lam22+2.0*lam01*lam02*lam12-lam11*lam02*lam02) + (lam12*lam12-lam11*lam22+lamB11*lamB22)*(lam22*(lam00*lamB11+lam11*lamB00)+lamB22*(lam00*lam11-lamB00*lamB11)-lamB00*lam12*lam12-lamB22*lam01*lam01-lamB11*lam02*lam02));
    return numerator/denominator;

}
//==============================================================================

double Transport::parallelThermalDiffusionCoefficient()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    // Get thermodynamic properties
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double me  = m_thermo.speciesMw(0)/NA;
    const double *const X = m_thermo.X();

    // Get collision integral information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar     = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute the lambdas
    double fac;
    double fac2 = 0.0;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;

    for (int i = 1; i < ns; ++i) {
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];

    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));

    // Second order
    return (lam22/(lam11*lam22-lam12*lam12));
}

//==============================================================================

double Transport::perpThermalDiffusionCoefficient()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    // Get thermodynamic properties
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double me  = m_thermo.speciesMw(0)/NA;
    const double *const X = m_thermo.X();
    const double P   = m_thermo.P();
    const double B = m_thermo.getBField();

    // Get collision integral information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute the lambdas
    double fac;
    double fac2 = 0.0;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;

    double lamB00;
    double lamB11;
    double lamB22;

    for (int i = 1; i < ns; ++i) {
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));

    lamB00 = QE*B/(KB*Te*fac2);
    lamB11 = 2.5*lamB00;
    lamB22 = 1.75*lamB11;

    // Second Order
    double numerator = lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11);
    double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);

    return numerator/denominator;
}

//==============================================================================

double Transport::transverseThermalDiffusionCoefficient()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    // Get thermodynamic properties
    const int ns     = m_thermo.nSpecies();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();
    const double nd  = m_thermo.numberDensity();
    const double me  = m_thermo.speciesMw(0)/NA;
    const double *const X = m_thermo.X();
    const double P   = m_thermo.P();
    const double B = m_thermo.getBField();

    // Get collision integral information
    const RealSymMat& Q11   = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22   = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute the lambdas
    double fac;
    double fac2 = 0.0;
    double lam11 = 0.0;
    double lam12 = 0.0;
    double lam22 = 0.0;

    double lamB00;
    double lamB11;
    double lamB22;

    for (int i = 1; i < ns; ++i) {
        lam11 += X[i]*Q11(i)*(25.0/4.0-3.0*Bstar(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));

    lamB00 = QE*B/(KB*Te*fac2);
    lamB11 = 2.5*lamB00;
    lamB22 = 1.75*lamB11;

    // Second Order
    double numerator = lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11);
    double denominator = (lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11);

    return numerator/denominator;
}
//==============================================================================

double Transport::sigmaParallel()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double* const X = m_thermo.X();
    const double nd = m_thermo.numberDensity();
    const double P = m_thermo.P();


    return 4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te)*(nd*KB*Te/P)* parallelDiffusionCoefficient();

}
//==============================================================================

double Transport::sigmaPerpendicular()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double* const X = m_thermo.X();
    const double nd = m_thermo.numberDensity();
    const double P = m_thermo.P();

    return 4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te)*(nd*KB*Te/P)* perpDiffusionCoefficient();

}

//==============================================================================

double Transport::sigmaTransverse()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double* const X = m_thermo.X();
    const double nd = m_thermo.numberDensity();
    const double P = m_thermo.P();

    return -4.0/25.0*(X[0]*X[0]*QE*QE)/(KB*KB*Te)*(nd*KB*Te/P)* transverseDiffusionCoefficient();

}

//==============================================================================

double Transport::parallelElectronThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    const double* const X = m_thermo.X();
    const double nd  = m_thermo.numberDensity();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();

    double Kt_parallel = parallelThermalDiffusionCoefficient();
    return (X[0]*X[0])*Kt_parallel;
}

//==============================================================================

double Transport::perpElectronThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

	if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    const double* const X = m_thermo.X();
    const double nd  = m_thermo.numberDensity();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();


    return (X[0]*X[0])*perpThermalDiffusionCoefficient();
}

//==============================================================================

double Transport::transverseElectronThermalConductivity()
{
	ERROR_IF_INTEGRALS_ARE_NOT_LOADED(0.0)

    if (!m_thermo.hasElectrons() || m_thermo.X()[0] < 1.0e-30)
        return 0.0;

    const double* const X = m_thermo.X();
    const double nd  = m_thermo.numberDensity();
    const double Th  = m_thermo.T();
    const double Te  = m_thermo.Te();

    return -(X[0]*X[0])*transverseThermalDiffusionCoefficient();
}

//==============================================================================

double Transport::ratioSigmaPerpPar()
{
    return sigmaPerpendicular()/sigmaParallel();
}
//==============================================================================

double Transport::ratioSigmaTransPar()
{
    return sigmaTransverse()/sigmaParallel();
}
//==============================================================================

double Transport::ratioLambdaPerpPar()
{
    return perpElectronThermalConductivity()/parallelElectronThermalConductivity();
}
//==============================================================================

double Transport::ratioLambdaTransPar()
{
    return transverseElectronThermalConductivity()/parallelElectronThermalConductivity();
}
//==============================================================================

std::vector<double> Transport::parallelThermalDiffusionRatio()
{
    double* species_values = new double [m_thermo.nSpecies()];

    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;

    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar   = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei = mp_collisions->Cstei(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute lambdas
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = 0.0;
    double lam02 = 0.0;
    std::vector<double> lam01ei(ns);
    std::vector<double> lam02ei(ns);
    double lam12 = 0.0;
    double lam22 = 0.0;
    double fac   = 0.0;


    for (int i = 1; i < ns; ++i) {
        fac = X[i]*Q11(i);
        lam00 += fac;
        lam01 += fac*(2.5-3.0*Cei(i));
        lam11 += fac*(25.0/4.0-3.0*Bstar(i));
        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }
    for (int i = 0; i < ns; ++i) {
    lam01ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(2.5*Q11(i)-3.0*Q12ei(i));
    lam02ei[i] =64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(35.0/8.0*Q11(i) - 21.0/2.0*Q12ei(i) + 6.0*Q13ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam02 = fac*lam02;
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));

    std::vector<double> kTi(ns);
    // Second Order
   /* kTi[0] = 2.5*Te/Th*X[0]*lam01/lam11;
    for (int i = 1; i < ns; ++i){
        kTi[i] = -2.5*Te/Th*X[0]*lam01ei[i]/lam11;
    }
    */

    // Third Order
    kTi[0] = 2.5*Te/Th*X[0]*(lam01*lam22 - lam02*lam12)/(lam11*lam22 - lam12*lam12);
    for (int i = 1; i < ns; ++i){
        kTi[i] = -2.5*Te/Th*X[0]*(lam01ei[i]*lam22 - lam02ei[i]*lam12)/(lam11*lam22 - lam12*lam12);
    }
    return kTi;

}
//==============================================================================
std::vector<double> Transport::perpThermalDiffusionRatio()
{
    // Get thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;
    const double B = m_thermo.getBField();

    // Get collision integrals
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar   = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei = mp_collisions->Cstei(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute lambdas
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = 0.0;
    double lam02 = 0.0;
    std::vector<double> lam01ei(ns);
    std::vector<double> lam02ei(ns);
    double lam12 = 0.0;
    double lam22 = 0.0;
    double fac   = 0.0;
    double fac2 = 0.0;

    double lamB00 = 0.0;
    double lamB11 = 0.0;
    double lamB22 = 0.0;

    for (int i = 1; i < ns; ++i) {
        fac = X[i]*Q11(i);
        lam00 += fac;
        lam01 += fac*(2.5-3.0*Cei(i));
        lam11 += fac*(25.0/4.0-3.0*Bstar(i));
        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));

    for (int i = 0; i < ns; ++i) {
        lam01ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(2.5*Q11(i)-3.0*Q12ei(i));
        lam02ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(35.0/8.0*Q11(i) - 21.0/2.0*Q12ei(i) + 6.0*Q13ei(i));
    }

    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam02 = fac*lam02;
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));

    lamB00 = QE*B/(KB*Te*fac2);
    lamB11 = 2.5*lamB00;
    lamB22 = 1.75*lamB11;

    std::vector<double> kTi(ns);

    // Third Order

    kTi[0] = 2.5*Te/Th*X[0]*(lam01*(lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11)) - lam02*lam12*(lam11*lam22-lamB11*lamB22-lam12*lam12)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));

    for (int i = 1; i < ns; ++i){
        kTi[i] = -2.5*Te/Th*X[0]*(lam01ei[i]*(lam22*(lam11*lam22-lamB11*lamB22-lam12*lam12) + lamB22*(lam11*lamB22+lam22*lamB11)) - lam02ei[i]*lam12*(lam11*lam22-lamB11*lamB22-lam12*lam12)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
    }
    return kTi;

}
//==============================================================================

std::vector<double> Transport::transverseThermalDiffusionRatio()
{
    // Get thermo properties
    const int ns = m_thermo.nSpecies();
    const double Th = m_thermo.T();
    const double Te = m_thermo.Te();
    const double nd = m_thermo.numberDensity();
    const double* const X = m_thermo.X();
    const double me = m_thermo.speciesMw(0)/NA;
    const double B = m_thermo.getBField();

    // Get collision integrals
    const RealSymMat& Q11 = mp_collisions->Q11(Th, Te, nd, X);
    const RealSymMat& Q22 = mp_collisions->Q22(Th, Te, nd, X);
    const RealSymMat& Bstar   = mp_collisions->Bstar(Th, Te, nd, X);
    const RealVector& Cei = mp_collisions->Cstei(Th, Te, nd, X);
    const RealVector& Q12ei = mp_collisions->Q12ei(Th, Te, nd, X);
    const RealVector& Q13ei = mp_collisions->Q13ei(Th, Te, nd, X);
    const RealVector& Q14ei = mp_collisions->Q14ei(Th, Te, nd, X);
    const RealVector& Q15ei = mp_collisions->Q15ei(Th, Te, nd, X);
    const double      Q23ee = mp_collisions->Q23ee(Th, Te, nd, X);
    const double      Q24ee = mp_collisions->Q24ee(Th, Te, nd, X);

    // Compute lambdas
    double lam00 = 0.0;
    double lam01 = 0.0;
    double lam11 = 0.0;
    double lam02 = 0.0;
    std::vector<double> lam01ei(ns);
    std::vector<double> lam02ei(ns);
    double lam12 = 0.0;
    double lam22 = 0.0;
    double fac   = 0.0;
    double fac2 = 0.0;

    double lamB00 = 0.0;
    double lamB11 = 0.0;
    double lamB22 = 0.0;

    for (int i = 1; i < ns; ++i) {
        fac = X[i]*Q11(i);
        lam00 += fac;
        lam01 += fac*(2.5-3.0*Cei(i));
        lam11 += fac*(25.0/4.0-3.0*Bstar(i));
        lam02 += X[i]*(35.0/8.0*Q11(i)-21.0/2.0*Q12ei(i)+6*Q13ei(i));
        lam12 += X[i]*(175.0/16.0*Q11(i)-315.0/8.0*Q12ei(i)+57.0*Q13ei(i)-30.0*
                       Q14ei(i));
        lam22 += X[i]*(1225.0/64.0*Q11(i)-735.0/8.0*Q12ei(i)+399.0/2.0*Q13ei(i)-
                       210.0*Q14ei(i)+90.0*Q15ei(i));
    }

    fac = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0];
    fac2 = 25.0/4.0*nd*KB*(Th/(X[0]*Te)+(1.0-Th/Te));

    for (int i = 0; i < ns; ++i) {
        lam01ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(2.5*Q11(i)-3.0*Q12ei(i));
        lam02ei[i] = 64.0/75.0*std::sqrt(me/(TWOPI*KB*KB*KB*Te))*X[0]*X[i]*(35.0/8.0*Q11(i) - 21.0/2.0*Q12ei(i) + 6.0*Q13ei(i));
    }

    lam00 = fac*lam00;
    lam01 = fac*lam01;
    lam11 = fac*(lam11 + SQRT2*X[0]*Q22(0));
    lam02 = fac*lam02;
    lam12 = fac*(lam12 + SQRT2*X[0]*(7.0/4.0*Q22(0)-2.0*Q23ee));
    lam22 = fac*(lam22 + SQRT2*X[0]*(77.0/16.0*Q22(0)-7.0*Q23ee+5.0*Q24ee));

    lamB00 = QE*B/(KB*Te*fac2);
    lamB11 = 2.5*lamB00;
    lamB22 = 1.75*lamB11;

    std::vector<double> kTi(ns);

    // Third Order

    kTi[0] = -2.5*Te/Th*X[0]*(lam01*(lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11)) + lam02*lam12*(lam11*lamB22+lam22*lamB11)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));

    for (int i = 1; i < ns; ++i){
        kTi[i] = 2.5*Te/Th*X[0]*(lam01ei[i]*(lamB22*(lam11*lam22-lamB11*lamB22-lam12*lam12) - lam22*(lam11*lamB22+lam22*lamB11)) + lam02ei[i]*lam12*(lam11*lamB22+lam22*lamB11)) / ((lam11*lam22-lamB11*lamB22-lam12*lam12)*(lam11*lam22-lamB11*lamB22-lam12*lam12) + (lam11*lamB22+lam22*lamB11)*(lam11*lamB22+lam22*lamB11));
    }
    return kTi;

}
//==============================================================================
std::vector<double> Transport::ratiokTPerpPar()
{
    const int ns = m_thermo.nSpecies();

    std::vector<double> ratio(ns);
    for (int i = 0; i < ns; ++i){
        ratio[i] = perpThermalDiffusionRatio()[i]/parallelThermalDiffusionRatio()[i];
    }
    return ratio;
}
//==============================================================================
std::vector<double> Transport::ratiokTTransPar()
{
    const int ns = m_thermo.nSpecies();

    std::vector<double> ratio(ns);
    for (int i = 0; i < ns; ++i){
        ratio[i] = transverseThermalDiffusionRatio()[i]/parallelThermalDiffusionRatio()[i];
    }
    return ratio;
}
//==============================================================================

    } // namespace Transport
} // namespace Mutation

