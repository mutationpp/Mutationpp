/**
 * @file Thermodynamics.cpp
 *
 * @brief Implementation of the Thermodynamics class.
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

#include "Thermodynamics.h"
#include "StateModel.h"
#include "Utilities.h"
#include "Composition.h"
#include "SpeciesListDescriptor.h"
//#include "MultiPhaseEquilSolver.h"
//#include "ParticleRRHO.h"

#include <set>

using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

//==============================================================================

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

Thermodynamics::Thermodynamics(
    const string& species_descriptor,
    const string& thermo_db)
    : mp_work1(NULL), mp_work2(NULL), mp_wrkcp(NULL), mp_default_composition(NULL),
      m_has_electrons(false), m_natoms(0), m_nmolecules(0)
{
    // Load the thermodynamic database
    mp_thermodb = Config::Factory<ThermoDB>::create(thermo_db, 0);
    if (!mp_thermodb->load(species_descriptor)) {
        cout << "Did not load all required species... Exiting." << endl;
        exit(1);
    }

    // Store the species and element order information for easy access
    for (int i = 0; i < nElements(); ++i)
        m_element_indices[element(i).name()] = i;
    
    for (int i = 0; i < nSpecies(); ++i) {
        m_species_indices[species(i).name()] = i;
        m_natoms += (species(i).type() == ATOM ? 1 : 0);
        m_nmolecules += (species(i).type() == MOLECULE ? 1 : 0);
    }
    
    m_has_electrons = (species(0).type() == ELECTRON);
    
    m_ngas = 0;
    for (int i = nSpecies()-1; (i > -1 && species(i).phase() != GAS); --i)
        m_ngas++;
    m_ngas = nSpecies()-m_ngas;

    // Build the composition matrix
    m_element_matrix = RealMatrix(nSpecies(), nElements());
    
    for (int i = 0; i < nSpecies(); ++i)
        for (int j = 0; j < nElements(); ++j)
            m_element_matrix(i,j) = 
                species(i).nAtoms(element(j).name());
    
    // Store the species molecular weights for faster access
    m_species_mw = RealVector(nSpecies());
    for (int i = 0; i < nSpecies(); ++i)
        m_species_mw(i) = species(i).molecularWeight();
    
    // Default composition (every element has equal parts)
    mp_default_composition = new double [nElements()];
    std::fill(
        mp_default_composition, mp_default_composition + nElements(), 
        1.0 / nElements());
    
    // Allocate a new equilibrium solver
    mp_equil = new MultiPhaseEquilSolver(*this);

    // Allocate storage for the work array
    mp_work1 = new double [nSpecies()];
    mp_work2 = new double [nSpecies()];
    mp_wrkcp = new double [nSpecies()*2]; /// @todo Need to fix all issues with Cp and Cv definitions!!!
    mp_y     = new double [nSpecies()];
}

//==============================================================================

Thermodynamics::~Thermodynamics()
{
    delete [] mp_work1;
    delete [] mp_work2;
    delete [] mp_wrkcp;
    delete [] mp_default_composition;
    delete [] mp_y;
    
    delete mp_thermodb;
    delete mp_equil;
}

//==============================================================================

void Thermodynamics::setDefaultComposition(const Composition& c)
{
    c.getComposition(m_element_indices, mp_default_composition);
    if (c.type() == Composition::MASS)
        convert<YE_TO_XE>(mp_default_composition, mp_default_composition);
}

//==============================================================================

int Thermodynamics::nPhases() const {
    return mp_equil->nPhases();
}

int Thermodynamics::nEnergyEqns() const {
    return stateModel()->nEnergyEqns();
}

int Thermodynamics::nMassEqns() const {
    return stateModel()->nMassEqns();
}

//==============================================================================

bool Thermodynamics::speciesThermoValidAtT(const size_t i, const double T) const
{
    return mp_thermodb->speciesThermoValidAtT(i, T);
}

//==============================================================================

void Thermodynamics::setState(
    const double* const p_v1, const double* const p_v2, const int vars)
{
    stateModel()->setState(p_v1, p_v2, vars);
    convert<X_TO_Y>(X(), mp_y);
}

//==============================================================================

double Thermodynamics::T() const {
    return stateModel()->T();
}

//==============================================================================

double Thermodynamics::Tr() const {
    return stateModel()->Tr();
}

//==============================================================================

double Thermodynamics::Tv() const {
    return stateModel()->Tv();
}

//==============================================================================

double Thermodynamics::Te() const {
    return stateModel()->Te();
}

//==============================================================================

double Thermodynamics::Tel() const {
    return stateModel()->Tel();
}

//==============================================================================

void Thermodynamics::getTemperatures(double* const p_T) const {
    stateModel()->getTemperatures(p_T);
}

//==============================================================================

void Thermodynamics::getEnergiesMass(double* const p_e) const{
    stateModel()->getEnergiesMass(p_e);
} 

//==============================================================================

void Thermodynamics::getEnthalpiesMass(double* const p_h) const {
    stateModel()->getEnthalpiesMass(p_h);
}

//==============================================================================

void Thermodynamics::getCpsMass(double* const p_cp) const {
    stateModel()->getCpsMass(p_cp);
}

//==============================================================================

void Thermodynamics::getCvsMass(double* const p_cv) const {
    stateModel()->getCvsMass(p_cv);
}

//==============================================================================

void Thermodynamics::getTagModes(int* const p_tag) const {
    stateModel()->getTagModes(p_tag);
}

//==============================================================================

double Thermodynamics::P() const {
    return stateModel()->P();
}

//==============================================================================

const double* const Thermodynamics::X() const {
    return stateModel()->X();
}

//==============================================================================

const double* const Thermodynamics::Y() const {
    return mp_y;
}

//==============================================================================

double Thermodynamics::standardStateT() const {
    return mp_thermodb->standardTemperature();
}

//==============================================================================

double Thermodynamics::standardStateP() const {
    return mp_thermodb->standardPressure();
}

//==============================================================================

double Thermodynamics::mixtureMw() const 
{
    return dot(Numerics::asVector(stateModel()->X(), nSpecies()), m_species_mw);
}

//==============================================================================

std::pair<int, int> Thermodynamics::equilibriumComposition(
    double T, double P, const double* const p_Xe, double* const p_X,
    MoleFracDef mdf) const
{
    return mp_equil->equilibrate(T, P, p_Xe, p_X, mdf);
}

//==============================================================================

void Thermodynamics::addEquilibriumConstraint(const double* const p_A)
{
    mp_equil->addConstraint(p_A);
}
    
//==============================================================================

void Thermodynamics::clearEquilibriumContraints()
{
    mp_equil->clearConstraints();
}

//==============================================================================

void Thermodynamics::elementPotentials(double *const p_lambda)
{
    mp_equil->elementPotentials(p_lambda);
}

//==============================================================================

void Thermodynamics::phaseMoles(double *const p_moles)
{
    mp_equil->phaseMoles(p_moles);
}

//==============================================================================

int Thermodynamics::nEquilibriumSteps() const {
    return mp_equil->nSteps();
}

//==============================================================================

int Thermodynamics::nEquilibriumNewtons() const {
    return mp_equil->nNewtons();
}

//==============================================================================

double Thermodynamics::numberDensity(const double T, const double P) const
{
    return P / (KB * T);
}

//==============================================================================

double Thermodynamics::numberDensity() const 
{
    double Xe = (hasElectrons() ? stateModel()->X()[0] : 0.0);
    double Th = stateModel()->T();
    double Te = stateModel()->Te();
    double P  = stateModel()->P();
    return P / (KB * (Th + Xe*(Te - Th)));
}

//==============================================================================

double Thermodynamics::pressure(
    const double T, const double rho, const double *const Y) const
{
    double pressure = 0.0;
    for (int i = 0; i < nGas(); ++i)
        pressure += Y[i] / speciesMw(i);
    pressure *= rho * T * RU;
    return pressure;
}

//==============================================================================
    
double Thermodynamics::density(
    const double T, const double P, const double *const X) const
{
    double density = 0.0;
    for (int i = 0; i < nGas(); ++i)
        density += X[i] * speciesMw(i);
    density *= P / (RU * T);
    return density;
}

//==============================================================================

double Thermodynamics::density() const
{
    return numberDensity() * mixtureMw() / NA;
    //return density(stateModel()->T(), stateModel()->P(), stateModel()->X());
}

//==============================================================================

void Thermodynamics::speciesCpOverR(double *const p_cp) const
{
    // WARNING: The total cp only makes sense at thermal equilibirum
    mp_thermodb->cp(
        stateModel()->T(), stateModel()->T(), stateModel()->T(), stateModel()->T(),
        stateModel()->T(), p_cp, NULL, NULL, NULL, NULL);
}

//==============================================================================

void Thermodynamics::speciesCpOverR(double T, double* const p_cp) const
{
    mp_thermodb->cp(
        T, T, T, T, T, p_cp, NULL, NULL, NULL, NULL);
}

//==============================================================================

void Thermodynamics::speciesCpOverR(
    double Th, double Te, double Tr, double Tv, double Tel, double *const p_cp, 
    double *const p_cpt, double *const p_cpr, double *const p_cpv, 
    double *const p_cpel) const
{
    mp_thermodb->cp(
        Th, Te, Tr, Tv, Tel, p_cp, p_cpt, p_cpr, p_cpv, p_cpel);
}

//==============================================================================

void Thermodynamics::speciesCvOverR(
    double Th, double Te, double Tr, double Tv, double Tel, double *const p_cv, 
    double *const p_cvt, double *const p_cvr, double *const p_cvv, 
    double *const p_cvel) const
{
    mp_thermodb->cp(Th, Te, Tr, Tv, Tel, p_cv, p_cvt, p_cvr, p_cvv, p_cvel);

    // WARNING: the total cv only makes sense at thermal equilibrium    
    if (p_cv != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cv[i] -= 1.0;
    }
    
    if (p_cvt != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cvt[i] -= 1.0;
    }
    
}

//==============================================================================

double Thermodynamics::mixtureFrozenCpMole() const
{
    return mixtureFrozenCpMass() * mixtureMw();
//    double cp = 0.0;
//    speciesCpOverR(mp_work1);
//    for (int i = 0; i < nGas(); ++i)
//        cp += mp_work1[i] * X()[i];
//    return (cp * RU);
}

//==============================================================================

double Thermodynamics::mixtureFrozenCpMass() const 
{
    double cp = 0.0;
    getCpsMass(mp_wrkcp);
    for (int i = 0; i < nSpecies(); ++i)
       cp += mp_wrkcp[i] * Y()[i];
    return cp;
//    return mixtureFrozenCpMole() / mixtureMw();
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumCpMole()
{
    const double T = this->T();
    
    // Compute species enthalpies and dg/dT
    speciesHOverRT(mp_work1);
    for (int j = 0; j < nSpecies(); ++j)
        mp_work2[j] = -mp_work1[j] / T;

    // Compute change in species moles due to change in temperature
    mp_equil->dNdg(mp_work2, mp_work2);

    // Now compute the sums which require dN/dT
    double sum1 = 0.0, sum2 = 0.0;
    for (int j = 0; j < nSpecies(); ++j) {
        sum1 += mp_work2[j];               // sum_j dN_j/dT
        sum2 += mp_work2[j]*mp_work1[j];   // sum_j dN_j/dT*H_j/RT
    }

    // Compute the species moles vector from the equilibrium solver
    mp_equil->speciesMoles(mp_work2);

    // Compute the sums that require N and H
    double sum3 = 0.0, sum4 = 0.0;
    for (int j = 0; j < nSpecies(); ++j) {
        sum3 += mp_work2[j];               // sum_j N_j
        sum4 += mp_work2[j]*mp_work1[j];   // sum_j N_j*H_j/RT
    }

    // Compute the Cp vector
    speciesCpOverR(mp_work1);
    double sum5 = 0.0;
    for (int j = 0; j < nSpecies(); ++j)
        sum5 += mp_work2[j]*mp_work1[j];   // sum_j N_j*Cp_j/R

    // Put together all the terms
    return RU*(sum3*(sum5+T*sum2)-T*sum1*sum4)/(sum3*sum3);
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumCpMass()
{
    const double T = this->T();

    // Compute species enthalpies and dg/dT
    speciesHOverRT(mp_work1);
    for (int j = 0; j < nSpecies(); ++j)
        mp_work2[j] = -mp_work1[j] / T;

    // Compute change in species moles due to change in temperature
    mp_equil->dNdg(mp_work2, mp_work2);

    // Now compute the sums which require dN/dT
    double sum1 = 0.0, sum2 = 0.0;
    for (int j = 0; j < nSpecies(); ++j) {
        sum1 += mp_work2[j]*speciesMw(j);  // sum_j dN_j/dT*Mw_j
        sum2 += mp_work2[j]*mp_work1[j];   // sum_j dN_j/dT*H_j/RT
    }

    // Compute the species moles vector from the equilibrium solver
    mp_equil->speciesMoles(mp_work2);

    // Compute the sums that require N and H
    double sum3 = 0.0, sum4 = 0.0;
    for (int j = 0; j < nSpecies(); ++j) {
        sum3 += mp_work2[j]*speciesMw(j);  // sum_j N_j
        sum4 += mp_work2[j]*mp_work1[j];   // sum_j N_j*H_j/RT
    }

    // Compute the Cp vector
    speciesCpOverR(mp_work1);
    double sum5 = 0.0;
    for (int j = 0; j < nSpecies(); ++j)
        sum5 += mp_work2[j]*mp_work1[j];   // sum_j N_j*Cp_j/R

    // Put together all the terms
    return RU*(sum3*(sum5+T*sum2)-T*sum1*sum4)/(sum3*sum3);
}

//==============================================================================

void Thermodynamics::dXidT(double* const p_dxdt) const
{   
    const double T = this->T();

    // Compute species enthalpies and dg/dT
    speciesHOverRT(p_dxdt);
    for (int j = 0; j < nSpecies(); ++j)
        p_dxdt[j] = -p_dxdt[j] / T;

    mp_equil->dXdg(p_dxdt, p_dxdt);
}

void Thermodynamics::dXidP(double* const p_dxdp) const
{
    const double P = this->P();
    for (int i = 0; i < nGas(); ++i)
        p_dxdp[i] = 1.0/P;
    for (int i = nGas(); i < nSpecies(); ++i)
        p_dxdp[i] = 0.0;
    mp_equil->dXdg(p_dxdp, p_dxdp);
}

//==============================================================================

double Thermodynamics::dRhodP()
{    
    // Get rho, P, T
    const double rho = density();
    const double P = this->P();
    const double Mwmix = mixtureMw();

    // Compute dX/dP (work2)
    for (int i = 0; i < nGas(); ++i)
        mp_work2[i] = 1.0/P;
    for (int i = nGas(); i < nSpecies(); ++i)
        mp_work2[i] = 0.0;
    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dT
    double dMwdP = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdP += mp_work2[i] * speciesMw(i);

    return (rho*(dMwdP/Mwmix+1.0/P));
}

//==============================================================================

double Thermodynamics::mixtureFrozenCvMole() const
{
    return mixtureFrozenCpMole() - RU;
}

//==============================================================================

double Thermodynamics::mixtureFrozenCvMass() const
{
    return (mixtureFrozenCpMole() - RU) / mixtureMw();
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumCvMass()
{
    // Get rho, P, T
    const double rho = density();
    const double P = this->P();
    const double T = this->T();
    const double Mwmix = mixtureMw();
    const double* const p_X = X();

    // Compute h/RT
    speciesHOverRT(mp_work1);
    
    // Compute dX/dT
    for (int j = 0; j < nSpecies(); ++j)
        mp_work2[j] = -mp_work1[j] / T;
    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dT
    double dMwdT = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdT += mp_work2[i] * speciesMw(i);

    // Compute reactive Cp
    double cp = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work1[i] * (mp_work2[i]*Mwmix - p_X[i]*dMwdT);
    cp *= (T / Mwmix);

    // Add Frozen Cp
    speciesCpOverR(mp_work2);
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work2[i] * p_X[i];
    cp *= RU / Mwmix;

    // Compute dX/dP (work2)
    for (int i = 0; i < nGas(); ++i)
        mp_work2[i] = 1.0/P;
    for (int i = nGas(); i < nSpecies(); ++i)
        mp_work2[i] = 0.0;
    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dP
    double dMwdP = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdP += mp_work2[i] * speciesMw(i);

    // Compute de/dP
    double dedp = 0.0;
    for (int i = 0; i < nSpecies(); i++)
        dedp += (mp_work1[i] - 1.0) * (mp_work2[i]*Mwmix - p_X[i]*dMwdP);
    dedp *= (RU * T / (Mwmix*Mwmix));

    // Compute density and energy derivatives
    const double drdt = dMwdT/Mwmix - 1.0/T; // note we leave out rho*(...)
    const double drdp = dMwdP/Mwmix + 1.0/P; // here also

    return (cp + (P/rho - dedp/drdp)*drdt);
}

void Thermodynamics::mixtureEnergies(double* const p_e) const
{
    stateModel()->getMixtureEnergiesMass(p_e);
}


/*double Thermodynamics::mixtureEquilibriumCvMass()
{
    // We assume that Y is already in equilibrium at T and P
    // Need to perturb current equilibrium conditions and compute dh/dT, dh/dP,
    // drho/T and drho/P
    
    double rho, rhoP, rhoT, T, P;
    double e, eP, eT, dT, dP;
    double drhodt, drhodp, dedt, dedp;
    
    // Get current mixture properties
    rho = density();
    e = mixtureEnergyMass();
    T = this->T();
    P = this->P();
    elementFractions(X(), mp_work2);
    
    // Perturb pressure
    dP = 1.0e-3*P;//std::sqrt(NumConst<double>::eps)*P;
    equilibrate(T, P+dP, mp_work2, mp_work1);
    rhoP = density();
    eP   = mixtureEnergyMass();
    
    // Perturb temperature
    dT = std::sqrt(NumConst<double>::eps)*T;
    equilibrate(T+dT, P, mp_work2, mp_work1);
    rhoT  = density();
    eT    = mixtureEnergyMass();
    
    // Return to current state
    equilibrate(T, P, mp_work2, mp_work1);
    drhodt = (rhoT-rho)/dT;
    drhodp = (rhoP-rho)/dP;
    dedt = (eT-e)/dT;
    dedp = (eP-e)/dP;
    
    //cout << endl << "de/dT: " << dedt << endl << "de/dP: " << dedp << endl;
    //cout << "dr/dT: " << drhodt << endl << "dr/dP: " << drhodp << endl;
    
    //cout << endl;
    //cout << T << " " << dT << " " << P << " " << dP << endl;
    //cout << drhodt << " " << drhodp << " " << dedt << " " << dedp << endl;
    
    // Compute Cv
    return (dedt - dedp*drhodt/drhodp);
}*/

//==============================================================================
    
double Thermodynamics::mixtureFrozenGamma() const 
{
    double cp = mixtureFrozenCpMole();
    return (cp / (cp - RU));
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumGamma()
{
    // Get rho, P, T
    const double rho = density();
    const double P = this->P();
    const double T = this->T();
    const double Mwmix = mixtureMw();
    const double* const p_X = X();

    speciesHOverRT(mp_work1);
    for (int j = 0; j < nSpecies(); ++j)
        mp_work2[j] = -mp_work1[j] / T;

    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dT
    double dMwdT = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdT += mp_work2[i] * speciesMw(i);

    // Compute reactive Cp
    double cp = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work1[i] * (mp_work2[i]*Mwmix - p_X[i]*dMwdT);
    cp *= (T / Mwmix);

    // Add Frozen Cp
    speciesCpOverR(mp_work2);
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work2[i] * p_X[i];
    cp *= RU / Mwmix;

    // Compute dX/dP (work2)
    for (int i = 0; i < nGas(); ++i)
        mp_work2[i] = 1.0/P;
    for (int i = nGas(); i < nSpecies(); ++i)
        mp_work2[i] = 0.0;
    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dP
    double dMwdP = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdP += mp_work2[i] * speciesMw(i);

    // Compute de/dP
    double dedp = 0.0;
    for (int i = 0; i < nSpecies(); i++)
        dedp += (mp_work1[i] - 1.0) * (mp_work2[i]*Mwmix - p_X[i]*dMwdP);
    dedp *= (RU * T / (Mwmix*Mwmix));

    // Compute density and energy derivatives
    const double drdt = dMwdT/Mwmix - 1.0/T; // note we leave out rho*(...)
    const double drdp = dMwdP/Mwmix + 1.0/P; // here also

    //return (dedt - dedp*drdt/drdp);
    return cp / (cp + (P/rho - dedp/drdp)*drdt);
}

//==============================================================================

double Thermodynamics::equilibriumSoundSpeed()
{
    // Get rho, P, T
    const double rho = density();
    const double P = this->P();
    const double T = this->T();
    const double Mwmix = mixtureMw();
    const double* const p_X = X();

    speciesHOverRT(mp_work1);
    for (int j = 0; j < nSpecies(); ++j)
        mp_work2[j] = -mp_work1[j] / T;

    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dT
    double dMwdT = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdT += mp_work2[i] * speciesMw(i);

    // Compute reactive Cp
    double cp = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work1[i] * (mp_work2[i]*Mwmix - p_X[i]*dMwdT);
    cp *= (T / Mwmix);

    // Add Frozen Cp
    speciesCpOverR(mp_work2);
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work2[i] * p_X[i];
    cp *= RU / Mwmix;

    // Compute dX/dP (work2)
    for (int i = 0; i < nGas(); ++i)
        mp_work2[i] = 1.0/P;
    for (int i = nGas(); i < nSpecies(); ++i)
        mp_work2[i] = 0.0;
    mp_equil->dXdg(mp_work2, mp_work2);

    // Compute dMw/dP
    double dMwdP = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        dMwdP += mp_work2[i] * speciesMw(i);

    // Compute de/dP
    double dedp = 0.0;
    for (int i = 0; i < nSpecies(); i++)
        dedp += (mp_work1[i] - 1.0) * (mp_work2[i]*Mwmix - p_X[i]*dMwdP);
    dedp *= (RU * T / (Mwmix*Mwmix));

    // Compute density and energy derivatives
    const double drdt = dMwdT/Mwmix - 1.0/T; // note we leave out rho*(...)
    const double drdp = dMwdP/Mwmix + 1.0/P; // here also

    return std::sqrt(cp / ((cp + (P/rho - dedp/drdp)*drdt) * rho * drdp));
}

//==============================================================================

void Thermodynamics::speciesHOverRT(
    double* const h, double* const ht, double* const hr, double* const hv,
    double* const hel, double* const hf) const
{
    mp_thermodb->enthalpy(
        stateModel()->T(), stateModel()->Te(), stateModel()->Tr(), stateModel()->Tv(),
        stateModel()->Tel(), h, ht, hr, hv, hel, hf);
}

//==============================================================================

void Thermodynamics::speciesHOverRT(double T, double* const h) const 
{
    mp_thermodb->enthalpy(
        T, T, T, T ,T, h, NULL, NULL, NULL, NULL, NULL);
}

//==============================================================================
void Thermodynamics::speciesHOverRT(
    double T, double Te, double Tr, double Tv, double Tel,
    double* const h, double* const ht, double* const hr, double* const hv,
    double* const hel, double* const hf) const
{
    mp_thermodb->enthalpy(
        T, Te, Tr, Tv, Tel,
        h, ht, hr, hv, hel, hf);
}

//==============================================================================

double Thermodynamics::mixtureHMole() const
{
    double h = 0.0;
    speciesHOverRT(mp_work1);
    for (int i = 0; i < nSpecies(); ++i)
        h += mp_work1[i] * X()[i];
    return (h * RU * stateModel()->T());
}

//==============================================================================

double Thermodynamics::mixtureHMass() const 
{
    return mixtureHMole() / mixtureMw();
}

//==============================================================================

void Thermodynamics::speciesSOverR(double *const p_s) const
{
    mp_thermodb->entropy(
        stateModel()->T(), stateModel()->Te(), stateModel()->Tr(), stateModel()->Tv(),
        stateModel()->Tel(), standardStateP(), p_s, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(stateModel()->P() / standardStateP());
    for (int i = 0; i < nSpecies(); ++i)
        if (species(i).phase() == GAS)
            p_s[i] -= lnp;
}

//==============================================================================

double Thermodynamics::mixtureSMole() const 
{
    double s = 0;
    speciesSOverR(mp_work1);
    for (int i = 0; i < nSpecies(); ++i)
        s += (mp_work1[i] - std::log(X()[i])) * X()[i];
    return (s * RU);
}

//==============================================================================

double Thermodynamics::mixtureSMass() const {
    return mixtureSMole() / mixtureMw();
}

//==============================================================================

void Thermodynamics::speciesGOverRT(double* const p_g) const
{
    mp_thermodb->gibbs(
        stateModel()->T(), stateModel()->Te(), stateModel()->Tr(), stateModel()->Tv(),
        stateModel()->Tel(), standardStateP(), p_g, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(stateModel()->P() / standardStateP());
    for (int i = 0; i < nSpecies(); ++i)
        if (species(i).phase() == GAS)
            p_g[i] += lnp;
}

//==============================================================================

void Thermodynamics::speciesGOverRT(double T, double P, double* const p_g) const
{
    mp_thermodb->gibbs(T, T, T, T, T, standardStateP(), p_g, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(P / standardStateP());
    for (int i = 0; i < nSpecies(); ++i)
        if (species(i).phase() == GAS)
            p_g[i] += lnp;
}

//==============================================================================

void Thermodynamics::speciesSTGOverRT(double T, double* const p_g) const {
    mp_thermodb->gibbs(
        T, T, T, T, T, standardStateP(), p_g, NULL, NULL, NULL, NULL);
}

//==============================================================================

void Thermodynamics::elementMoles(
    const double *const species_N, double *const element_N) const
{
    asVector(element_N, nElements()) = 
        asVector(species_N, nSpecies()) * m_element_matrix;
}

//==============================================================================

void Thermodynamics::elementFractions(
    const double* const Xs, double* const Xe) const
{
//    RealVecWrapper wrapper(Xe, nElements());
//    wrapper = asVector(Xs, nSpecies()) * m_element_matrix;
//    double sum = wrapper.sum();
//    for (int i = 0; i < nElements(); ++i)
//        Xe[i] /= sum;
    //wrapper = wrapper / wrapper.sum();

    const int ne = nElements();
    const int ns = nSpecies();

    double temp = Xs[0];
    for (int k = 0; k < ne; ++k)
        Xe[k] = m_element_matrix(0,k)*temp;

    for (int i = 1; i < ns; ++i) {
        temp = Xs[i];
        for (int k = 0; k < ne; ++k)
            Xe[k] += m_element_matrix(i,k)*temp;
    }

    temp = 0.0;
    for (int k = 0; k < ne; ++k) {
        Xe[k] = std::max(0.0, Xe[k]);
        temp += Xe[k];
    }

    for (int k = 0; k < ne; ++k)
        Xe[k] /= temp;
}

//==============================================================================

void Thermodynamics::surfaceMassBalance(
    const double *const p_Yke, const double *const p_Ykg, const double T, 
    const double P, const double Bg, double &Bc, double &hw, double *const p_Xs)
{
    const int ne = nElements();
    const int ng = nGas();
    
    double p_Xw [ne];
    double* p_X  = (p_Xs != NULL ? p_Xs : mp_work1);
    double* p_h  = mp_work2;
    
    // Initialize the wall element fractions to be the pyrolysis gas fractions
    double sum = 0.0;
    for (int i = 0; i < ne; ++i) {
        p_Xw[i] = p_Yke[i] + Bg*p_Ykg[i];
        sum += p_Xw[i];
    }
    
    // Use "large" amount of carbon to simulate infinite char
    int ic = elementIndex("C");
    //double carbon = std::min(1000.0, std::max(100.0,1000.0*Bg));
    double carbon = std::max(100.0*Bg, 200.0);
    p_Xw[ic] += carbon;
    sum += carbon;
    
    for (int i = 0; i < ne; ++i)
        p_Xw[i] /= sum;
    
    // Compute equilibrium
    convert<YE_TO_XE>(p_Xw, p_Xw);
    equilibriumComposition(T, P, p_Xw, p_X, IN_PHASE);
    
    // Compute the gas mass fractions at the wall
    double mwg = 0.0;
    double ywc = 0.0;
    
    for (int j = 0; j < ng; ++j) {
        mwg += speciesMw(j) * p_X[j];
        ywc += elementMatrix()(j,ic) * p_X[j];
        //for (int i = 0; i < ne; ++i)
        //    p_Yw[i] += elementMatrix()(j,i) * p_X[j];
    }
    
    //for (int i = 0; i < ne; ++i)
    //    p_Yw[i] *= atomicMass(i) / mwg;
    ywc *= atomicMass(ic) / mwg;
    
    // Compute char mass blowing rate
    Bc = (p_Yke[ic] + Bg*p_Ykg[ic] - ywc*(1.0 + Bg)) / (ywc - 1.0);
    Bc = std::max(Bc, 0.0);
    
    // Compute the gas enthalpy
    speciesHOverRT(T, p_h);
    
    hw = 0.0;
    for (int i = 0; i < ng; ++i)
        hw += p_X[i] * p_h[i];
    hw *= RU * T / mwg;
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation


