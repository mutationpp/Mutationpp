/**
 * @file Thermodynamics.cpp
 *
 * @brief Implementation of the Thermodynamics class.
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

#include "Thermodynamics.h"
#include "StateModel.h"
#include "Utilities.h"
#include "Composition.h"

#include <set>

using namespace std;
using namespace Mutation::Utilities;
using namespace Eigen;

//==============================================================================

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

Thermodynamics::Thermodynamics(
    const string& species_descriptor,
    const string& thermo_db,
    const string& state_model )
    : mp_work1(NULL), mp_work2(NULL), mp_wrkcp(NULL), mp_default_composition(NULL),
      m_has_electrons(false), m_natoms(0), m_nmolecules(0)
{
    try {
        // Load the thermodynamic database
        mp_thermodb = Config::Factory<ThermoDB>::create(thermo_db, 0);
    } catch (Error& e) {
        e << "\nWas trying to load the thermodynamic database.";
        throw;
    }
    if (!mp_thermodb->load(species_descriptor)) {
        throw InvalidInputError("species list", species_descriptor)
            << "Could not find all required species in the thermodynamic "
            << "database.";
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
    m_element_matrix.resize(nSpecies(), nElements());
    
    for (int i = 0; i < nSpecies(); ++i)
        for (int j = 0; j < nElements(); ++j)
            m_element_matrix(i,j) = 
                species(i).nAtoms(element(j).name());
    
    // Store the species molecular weights for faster access
    m_species_mw.resize(nSpecies());
    for (int i = 0; i < nSpecies(); ++i)
        m_species_mw(i) = species(i).molecularWeight();
    
    // Default composition (every element has equal parts)
    mp_default_composition = new double [nElements()];
    std::fill(
        mp_default_composition, mp_default_composition + nElements(), 
        1.0 / nElements());
    
    // Allocate a new equilibrium solver
    mp_equil = new MultiPhaseEquilSolver(*this);
    
    try {
        // Allocate a new state model
        mp_state = Config::Factory<StateModel>::create(state_model, *this);
    } catch (Error& e) {
        // Can't save this, fail gracefully
        delete mp_thermodb;
        delete mp_default_composition;
        delete mp_equil;

        // Add a bit more info to the error message
        e << "\nWas trying to load the requested state model.";
        throw;
    }

    // Allocate storage for the work array
    mp_work1 = new double [nSpecies()];
    mp_work2 = new double [nSpecies()];
    mp_wrkcp = new double [nSpecies()*nEnergyEqns()];
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
    delete mp_state;
}

//==============================================================================

void Thermodynamics::setDefaultComposition(const Composition& c)
{
    c.getComposition(m_element_indices, mp_default_composition);
    if (c.type() == Composition::MASS)
        convert<YE_TO_XE>(mp_default_composition, mp_default_composition);
}

void Thermodynamics::setDefaultComposition(const double* const p_vec)
{
    std::copy(p_vec, p_vec + nElements(), mp_default_composition);
}

//==============================================================================

int Thermodynamics::nPhases() const {
    return mp_equil->nPhases();
}

int Thermodynamics::nEnergyEqns() const {
    return mp_state->nEnergyEqns();
}

int Thermodynamics::nMassEqns() const {
    return mp_state->nMassEqns();
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
    mp_state->setState(p_v1, p_v2, vars);
    convert<X_TO_Y>(X(), mp_y);
}

//==============================================================================

void Thermodynamics::setBField(const double B) {
    mp_state->setBField(B);
}

//==============================================================================

double Thermodynamics::getBField() const {
    return mp_state->getBField();
}

//==============================================================================

double Thermodynamics::T() const {
    return mp_state->T();
}

//==============================================================================

double Thermodynamics::Tr() const {
    return mp_state->Tr();
}

//==============================================================================

double Thermodynamics::Tv() const {
    return mp_state->Tv();
}

//==============================================================================

double Thermodynamics::Te() const {
    return mp_state->Te();
}

//==============================================================================

double Thermodynamics::Tel() const {
    return mp_state->Tel();
}

//==============================================================================

void Thermodynamics::getTemperatures(double* const p_T) const {
    mp_state->getTemperatures(p_T);
}

//==============================================================================

void Thermodynamics::getEnergiesMass(double* const p_e) const{
    mp_state->getEnergiesMass(p_e);     
} 

//==============================================================================

void Thermodynamics::getEnthalpiesMass(double* const p_h) const {
    mp_state->getEnthalpiesMass(p_h);
}

//==============================================================================

void Thermodynamics::getCpsMass(double* const p_cp) const {
    mp_state->getCpsMass(p_cp);
}

//==============================================================================

void Thermodynamics::getCvsMass(double* const p_cv) const {
    mp_state->getCvsMass(p_cv);
}

//==============================================================================

void Thermodynamics::getTagModes(int* const p_tag) const {
    mp_state->getTagModes(p_tag);
}

//==============================================================================

double Thermodynamics::P() const {
    return mp_state->P();
}

//==============================================================================

const double* const Thermodynamics::X() const {
    return mp_state->X();
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
    //return dot(Numerics::asVector(mp_state->X(), nSpecies()), m_species_mw);
    return (Map<const ArrayXd>(mp_state->X(), nSpecies())*m_species_mw).sum();
}

//==============================================================================

std::pair<int, int> Thermodynamics::equilibriumComposition(
    double T, double P, const double* const p_Xe, double* const p_X,
    MoleFracDef mdf) const
{
    return mp_equil->equilibrate(T, P, p_Xe, p_X, mdf);
}

//==============================================================================

void Thermodynamics::equilibrate(double T, double P, double* const p_Xe) const
{
    mp_state->equilibrate(T, P, p_Xe);
    convert<X_TO_Y>(X(), mp_y);
}

//==============================================================================

void Thermodynamics::addEquilibriumConstraint(const double* const p_A)
{
    //mp_equil->addConstraint(p_A);
    throw NotImplementedError("Thermodynamics::addEquilibriumConstraint()");
}
    
//==============================================================================

void Thermodynamics::clearEquilibriumContraints()
{
    //mp_equil->clearConstraints();
    throw NotImplementedError("Thermodynamics::clearEquilibriumContraints()");
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
    double Xe = (hasElectrons() ? mp_state->X()[0] : 0.0);
    double Th = mp_state->T();
    double Te = mp_state->Te();
    double P  = mp_state->P();
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
    //return density(mp_state->T(), mp_state->P(), mp_state->X());
}

//==============================================================================

void Thermodynamics::densities(double* const p_rho) const
{
    double rho = density();
    for (int i = 0; i < nGas(); ++i)
        p_rho[i] = rho * Y()[i];
}

//==============================================================================

void Thermodynamics::speciesCpOverR(double *const p_cp) const
{
    // WARNING: The total cp only makes sense at thermal equilibirum
    mp_thermodb->cp(
        mp_state->T(), mp_state->T(), mp_state->T(), mp_state->T(),
        mp_state->T(), p_cp, NULL, NULL, NULL, NULL);
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
    if (nSpecies() == 1)
        return mixtureFrozenCpMole();
    
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
    if (nSpecies() == 1)
        return mixtureFrozenCpMass();

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

void Thermodynamics::dXjdci(int i, double* const p_dxdc) const
{
    mp_equil->dXdc(i, p_dxdc);
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
    return mixtureFrozenCvMass() * mixtureMw();
    //return mixtureFrozenCpMole() - RU; // Wrong in the T-Tv case !!!
}

//==============================================================================

double Thermodynamics::mixtureFrozenCvMass() const
{
    double cv = 0.0;
    getCvsMass(mp_wrkcp);
    for (int i = 0; i < nSpecies(); ++i)
       cv += mp_wrkcp[i] * Y()[i];
    return cv;
    //return (mixtureFrozenCpMole() - RU) / mixtureMw(); // Wrong in the T-Tv case !!!
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

//==============================================================================

void Thermodynamics::sumSpeciesMass(
    const double* const p_s, double* const p_e, const int size) const
{
    const int ns = nSpecies();

    Map<VectorXd>(p_e, size) =
        Map<const Matrix<double, -1, -1, RowMajor> >(p_s, size, ns) *
        Map<const VectorXd>(mp_y, ns);
}

//==============================================================================

void Thermodynamics::mixtureEnergies(double* const p_e)
{
    // Compute the species energies
    mp_state->getEnergiesMass(mp_wrkcp);
    sumSpeciesMass(mp_wrkcp, p_e, nEnergyEqns());
}

//==============================================================================


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
    double cv = mixtureFrozenCvMole();
    return (cp / cv);
    //return (cp / (cp - RU));
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
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), h, ht, hr, hv, hel, hf);
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
    return (h * RU * mp_state->T());
}

//==============================================================================

double Thermodynamics::mixtureHMole(double T) const
{
    double h = 0.0;
    speciesHOverRT(T, mp_work1);
    for (int i = 0; i < nSpecies(); ++i)
        h += mp_work1[i] * X()[i];
    return (h * RU * T);
}

//==============================================================================

double Thermodynamics::mixtureHMass() const 
{
    return mixtureHMole() / mixtureMw();
}

//==============================================================================

double Thermodynamics::mixtureHMass(double T) const
{
    return mixtureHMole(T) / mixtureMw();
}

//==============================================================================

void Thermodynamics::speciesSOverR(double *const p_s) const
{
    mp_thermodb->entropy(
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), standardStateP(), p_s, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(mp_state->P() / standardStateP());
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
        s += (X()[i] > 0.0 ? (mp_work1[i] - std::log(X()[i])) * X()[i] : 0.0);
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
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), standardStateP(), p_g, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(mp_state->P() / standardStateP());
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
    Map<VectorXd>(element_N, nElements()) =
        Map<const VectorXd>(species_N, nSpecies()) * m_element_matrix;
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


