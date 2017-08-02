/**
 * @file Kinetics.cpp
 *
 * @brief Implementation of Kinetics class.
 */

/*
 * Copyright 2014-2017 von Karman Institute for Fluid Dynamics (VKI)
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

#include "Kinetics.h"
#include "Constants.h"
#include "Utilities.h"

using namespace std;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;


namespace Mutation {
    namespace Kinetics {

using Mutation::Thermodynamics::Thermodynamics;

//==============================================================================

Kinetics::Kinetics(
    const Thermodynamics& thermo, string mechanism)
    : m_thermo(thermo),
      mp_rates(NULL),
      m_thirdbodies(thermo.nSpecies(), m_thermo.hasElectrons()),
      m_jacobian(thermo),
      mp_ropf(NULL),
      mp_ropb(NULL),
      mp_rop(NULL),
      mp_wdot(NULL)
{
    if (mechanism == "none")
        return;
    
    // Get the path to the mechanism file
    mechanism = databaseFileName(mechanism, "mechanisms");

    // Open the mechanism file as an XML document
    IO::XmlDocument doc(mechanism);        
    IO::XmlElement root = doc.root();
    
    if (root.tag() != "mechanism") {
        cout << "Root element in mechanism file " << mechanism
             << " is not of 'mechanism' type!";
        exit(1); 
    }
    
    // Now loop over all of the reaction nodes and add each reaction to the
    // corresponding data structure pieces
    IO::XmlElement::const_iterator iter = root.begin();
    for ( ; iter != root.end(); ++iter) {        
        if (iter->tag() == "reaction")
            addReaction(Reaction(*iter, thermo));
        else if (iter->tag() == "arrhenius_units")
            Arrhenius::setUnits(*iter);
    }
    
    // Setup the rate manager
    mp_rates = new RateManager(thermo.nSpecies(), m_reactions);
    
    // Finally close the reaction mechanism
    closeReactions(true);
}

Kinetics::~Kinetics()
{
    if (mp_rates != NULL)
        delete mp_rates;
    if (mp_ropf != NULL)
        delete [] mp_ropf;
    if (mp_ropb != NULL)
        delete [] mp_ropb;
    if (mp_rop != NULL)
        delete [] mp_rop;
    if (mp_wdot != NULL)
        delete [] mp_wdot;
}

//==============================================================================

void Kinetics::addReaction(const Reaction& reaction)
{
    // Add reaction to reaction list
    m_reactions.push_back(reaction);
    
    // Insert the reactants
    m_reactants.addReaction(nReactions()-1, reaction.reactants());
    
    // Insert products
    if (reaction.isReversible())
        m_rev_prods.addReaction(nReactions()-1, reaction.products());
    else
        m_irr_prods.addReaction(nReactions()-1, reaction.products());
    
    // Add thirdbodies if necessary
    if (reaction.isThirdbody())
        m_thirdbodies.addReaction(nReactions()-1, reaction.efficiencies());
    
    // Add the reaction to the jacobian managaer
    m_jacobian.addReaction(reaction);
}

//==============================================================================

void Kinetics::closeReactions(const bool validate_mechanism) 
{
    const size_t ns = m_thermo.nSpecies();
    
    // Validate the mechanism
    if (validate_mechanism) {
        bool is_valid = true;
        //cout << "Validating reaction mechanism..." << endl;
        
        // Check for duplicate reactions
        for (size_t i = 0; i < nReactions()-1; ++i) {
            for (size_t j = i+1; j < nReactions(); ++j) {
                if (m_reactions[i] == m_reactions[j]) {
                    cerr << "Reactions " << i+1 << " \"" 
                         << m_reactions[i].formula()
                         << "\" and " << j+1 << " \""
                         << m_reactions[j].formula()
                         << "\" are identical." << endl;
                    is_valid = false;
                }
            }
        }
        
        // Check for elemental mass and charge conservation
        for (size_t i = 0; i < nReactions(); ++i) {
            if (!m_reactions[i].conservesChargeAndMass()) {
                cerr << "Reaction " << i+1 << " \"" << m_reactions[i].formula()
                     << "\" does not conserve charge or mass." << endl;
                is_valid = false;
            }
        }
        
        if (!is_valid) {
            cout << "Validation checks failed!" << endl;
            exit(1);
        }
    }
    
    // Allocate work arrays
    mp_ropf  = new double [nReactions()];
    mp_ropb  = new double [nReactions()];
    mp_rop   = new double [nReactions()];
    mp_wdot  = new double [m_thermo.nSpecies()];
    
}

//==============================================================================

void Kinetics::getReactionDelta(
    const double* const p_s, double* const p_r) const
{
    if (nReactions() == 0)
        return;

    m_reactants.decrReactions(p_s, p_r);
    m_rev_prods.incrReactions(p_s, p_r);
    m_irr_prods.incrReactions(p_s, p_r);
}

//==============================================================================

/*vector<size_t> Kinetics::speciesIndices(
    const multiset<string>& set)
{    
    multiset<string>::const_iterator iter = set.begin();
    vector<size_t> indices;
    
    for ( ; iter != set.end(); ++iter)
        indices.push_back(m_thermo.speciesIndex(*iter));
    
    return indices;
}*/

//==============================================================================

/*vector<pair<size_t, double> > Kinetics::thirdbodyEffs(
    const vector<pair<string, double> >& string_effs)
{
    vector<pair<size_t, double> > effs;
    vector<pair<string, double> >::const_iterator iter;
    
    for (iter = string_effs.begin(); iter != string_effs.end(); ++iter)
        effs.push_back(
            make_pair(m_thermo.speciesIndex(iter->first), iter->second));
    
    return effs;
}*/

//==============================================================================

void Kinetics::forwardRateCoefficients(double* const p_kf)
{
    if (nReactions() == 0)
        return;

    mp_rates->update(m_thermo);
    const double* const p_lnkf = mp_rates->lnkf();
    for (int i = 0; i < nReactions(); ++i)
        p_kf[i] = std::exp(p_lnkf[i]);
}

//==============================================================================

void Kinetics::backwardRateCoefficients(double* const p_kb)
{
    if (nReactions() == 0)
        return;

    mp_rates->update(m_thermo);
    const double* const p_lnkb = mp_rates->lnkb();
    for (int i = 0; i < nReactions(); ++i)
        p_kb[i] = std::exp(p_lnkb[i]);
}

/*
//==============================================================================

void Kinetics::forwardRatesOfProgress(
    const double T, const double* const p_conc, double* const p_ropf)
{
    forwardRateCoefficients(T, p_ropf);
    m_reactants.multReactions(p_conc, p_ropf);
    m_thirdbodies.multiplyThirdbodies(p_conc, p_ropf);
}

//==============================================================================

void Kinetics::backwardRatesOfProgress(
    const double T, const double* const p_conc, double* const p_ropb)
{
    backwardRateCoefficients(T, p_ropb);
    m_rev_prods.multReactions(p_conc, p_ropb);
    m_thirdbodies.multiplyThirdbodies(p_conc, p_ropb);
}

//==============================================================================

void Kinetics::updateROP(
    const double T, const double* const p_conc, double* const p_rop)
{
    forwardRateCoefficients(T, mp_ropf);
    m_reactants.multReactions(p_conc, mp_ropf);
    backwardRateCoefficients(T, mp_ropb);
    m_rev_prods.multReactions(p_conc, mp_ropb);
    
    for (int i = 0; i < m_num_rxns; ++i)
        p_rop[i] = (mp_ropf[i] - mp_ropb[i]);
    
    m_thirdbodies.multiplyThirdbodies(p_conc, p_rop);
}*/

//==============================================================================

void Kinetics::netRatesOfProgress(double* const p_rop)
{
    // Special case of no reactions
    if (nReactions() == 0) {
        for (int i = 0; i < nReactions(); ++i)
            p_rop[i] = 0.0;
        return;
    }

    // Compute species concentrations (mol/m^3)
    const double mix_conc = m_thermo.numberDensity() / NA;
    const double* const p_x = m_thermo.X();
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        mp_wdot[i] = p_x[i] * mix_conc;
    
    // Update the forward and backward rate coefficients
    mp_rates->update(m_thermo);
    
    // Compute forward ROP
    const double* const lnkf = mp_rates->lnkf();
    for (int i = 0; i < nReactions(); ++i)
        mp_ropf[i] = std::exp(lnkf[i]);
    m_reactants.multReactions(mp_wdot, mp_ropf);
    
    // Compute reverse ROP
    const double* const lnkb = mp_rates->lnkb();
    for (int i = 0; i < nReactions(); ++i)
        mp_ropb[i] = std::exp(lnkb[i]);
    m_rev_prods.multReactions(mp_wdot, mp_ropb);
    
    // Compute net ROP
    for (int i = 0; i < nReactions(); ++i)
        mp_rop[i] = mp_ropf[i] - mp_ropb[i];
    
    // Thirdbody efficiencies
    m_thirdbodies.multiplyThirdbodies(mp_wdot, mp_rop);
    
    // Multiply by species molecular weights
    for (int i = 0; i < nReactions(); ++i)
        p_rop[i] = mp_rop[i];
}

/*
//==============================================================================

void Kinetics::netProductionRates(
    const double T, const double* const p_conc, double* const p_wdot)
{
    std::fill(p_wdot, p_wdot+m_thermo.nSpecies(), 0.0);

    netRatesOfProgress(T, p_conc, mp_rop);
    m_reactants.decrSpecies(mp_rop, p_wdot);
    m_rev_prods.incrSpecies(mp_rop, p_wdot);
    m_irr_prods.incrSpecies(mp_rop, p_wdot);
    
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_wdot[i] *= m_thermo.speciesMw(i);
}*/

//==============================================================================

void Kinetics::netProductionRates(double* const p_wdot)
{
    // Special case of no reactions
    if (nReactions() == 0) {
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            p_wdot[i] = 0.0;
        return;
    }

    // Compute species concentrations (mol/m^3)
    const double mix_conc = m_thermo.numberDensity() / NA;
    const double* const p_x = m_thermo.X();
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_wdot[i] = p_x[i] * mix_conc;

    
    // Update the forward and backward rate coefficients
    mp_rates->update(m_thermo);
    
    // Compute forward ROP
    //cout << m_thermo.T() << ", " << m_thermo.Te() << endl;
    const double* const lnkf = mp_rates->lnkf();
    //cout << "ln forward rate" << endl;
    for (int i = 0; i < nReactions(); ++i) {
    //    cout << lnkf[i] << endl;
        mp_ropf[i] = std::exp(lnkf[i]);
    }
    m_reactants.multReactions(p_wdot, mp_ropf);
    
    // Compute reverse ROP
    const double* const lnkb = mp_rates->lnkb();
    //cout << "ln backward rate" << endl;
    for (int i = 0; i < nReactions(); ++i) {
    //    cout << lnkb[i] << endl;
        mp_ropb[i] = std::exp(lnkb[i]);
    }
    m_rev_prods.multReactions(p_wdot, mp_ropb);
    
    // Compute net ROP
    for (int i = 0; i < nReactions(); ++i)
        mp_rop[i] = mp_ropf[i] - mp_ropb[i];

    // Thirdbody efficiencies
    m_thirdbodies.multiplyThirdbodies(p_wdot, mp_rop);
    
    // Sum all contributions from every reaction
    std::fill(p_wdot, p_wdot+m_thermo.nSpecies(), 0.0);
    m_reactants.decrSpecies(mp_rop, p_wdot);
    m_rev_prods.incrSpecies(mp_rop, p_wdot);
    m_irr_prods.incrSpecies(mp_rop, p_wdot);


    // Multiply by species molecular weights
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_wdot[i] *= m_thermo.speciesMw(i);
}

//==============================================================================

void Kinetics::jacobianRho(double* const p_jac)
{
    // Special case of no reactions
    if (nReactions() == 0) {
        for (int i = 0; i < m_thermo.nSpecies()*m_thermo.nSpecies(); ++i)
            p_jac[i] = 0.0;
        return;
    }

    // Update reaction rate coefficients
    mp_rates->update(m_thermo);
    
    const double* const lnkf = mp_rates->lnkf();
    for (int i = 0; i < nReactions(); ++i)
        mp_ropf[i] = std::exp(lnkf[i]);
    
    const double* const lnkb = mp_rates->lnkb();
    for (int i = 0; i < nReactions(); ++i)
        mp_ropb[i] = std::exp(lnkb[i]);
    
    // Compute species concentrations (mol/m^3)
    const double mix_conc = m_thermo.numberDensity() / NA;
    const double* const p_x = m_thermo.X();
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        mp_rop[i] = p_x[i] * mix_conc;
    
    // Compute the Jacobian matrix
    m_jacobian.computeJacobian(mp_ropf, mp_ropb, mp_rop, p_jac);
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation

