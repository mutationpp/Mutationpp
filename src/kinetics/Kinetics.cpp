/**
 * @file Kinetics.cpp
 *
 * @brief Implementation of Kinetics class.
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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
#include "Numerics.h"

#include <algorithm>
#include <functional>
#include <cmath>
#include <math.h>

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
using namespace Mutation::Numerics;

namespace Mutation {
    namespace Kinetics {

using Mutation::Thermodynamics::Thermodynamics;

//==============================================================================

Kinetics::Kinetics(
    const Thermodynamics& thermo, string mechanism)
    : m_name("unnamed"),
      m_thermo(thermo),
      mp_rates(NULL),
      m_thirdbodies(thermo.nSpecies(), m_thermo.hasElectrons()),
      m_jacobian(thermo),
      mp_ropf(NULL),
      mp_ropb(NULL),
      mp_rop(NULL),
      mp_wdot(NULL),
      mp_r(NULL),
      mp_r_22(NULL),
      mp_r_23(NULL)
{
    if (mechanism == "none")
        return;
    
    // Get the path to the mechanism file
    mechanism = databaseFileName(mechanism, "mechanisms");

    // Open the mechanism file as an XML document
    IO::XmlDocument doc(mechanism);        
    IO::XmlElement root = doc.root();
    
    if (root.tag() != "mechanism") {
        throw FileParseError(doc.file(), root.line())
            << "Root element in mechanism file " << mechanism
            << " is not of 'mechanism' type!";
    }
    
    // Get the mechanism name
    root.getAttribute("name", m_name, m_name);

    // Now loop over all of the reaction nodes and add each reaction to the
    // corresponding data structure pieces
    IO::XmlElement::const_iterator iter = root.begin();
    for ( ; iter != root.end(); ++iter) {        
        if (iter->tag() == "reaction")
            addReaction(Reaction(*iter, thermo));
        else if (iter->tag() == "arrhenius_units")
            Arrhenius::setUnits(*iter);
	//else if (iter->tag() == "g2t_units")
	//    G2T::setUnits(*iter);
	//else if (iter->tag() == "rationalexp_units")
	//    RationalExp::setUnits(*iter);
	//else if (iter->tag() == "constRate_units")
	//    constRate::setUnits(*iter);
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
    if (mp_r != NULL)
        delete [] mp_r;
    if (mp_r_22 != NULL)
        delete [] mp_r_22;
    if (mp_r_23 != NULL)
        delete [] mp_r_23;
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
        // Check for duplicate reactions
        for (size_t i = 0; i < nReactions()-1; ++i)
            for (size_t j = i+1; j < nReactions(); ++j)
                if (m_reactions[i] == m_reactions[j])
                    throw InvalidInputError("mechanism", m_name)
                        << "Reactions " << i+1 << " \""
                        << m_reactions[i].formula()
                        << "\" and " << j+1 << " \""
                        << m_reactions[j].formula()
                        << "\" are identical.";

        // Check for elemental mass and charge conservation
        for (size_t i = 0; i < nReactions(); ++i)
            if (!m_reactions[i].conservesChargeAndMass())
                throw InvalidInputError("mechanism", m_name)
                    << "Reaction " << i+1 << " \"" << m_reactions[i].formula()
                    << "\" does not conserve charge or mass.";
    }
    
    // Allocate work arrays
    mp_ropf  = new double [nReactions()];
    mp_ropb  = new double [nReactions()];
    mp_rop   = new double [nReactions()];
    mp_wdot  = new double [m_thermo.nSpecies()];

    mp_r     = new double [m_thermo.nMolecules()];
    mp_r_22  = new double [m_thermo.nMolecules()];
    mp_r_23  = new double [m_thermo.nMolecules()];
    
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
    Map<ArrayXd>(p_kf, nReactions()) = 
        Map<const ArrayXd>(mp_rates->lnkf(), nReactions()).exp();
}

//==============================================================================

void Kinetics::backwardRateCoefficients(double* const p_kb)
{
    if (nReactions() == 0)
        return;

    mp_rates->update(m_thermo);
    Map<ArrayXd>(p_kb, nReactions()) = 
        Map<const ArrayXd>(mp_rates->lnkb(), nReactions()).exp();
        
    for(int i=0; i < mp_rates->irrReactions().size(); ++i)
        p_kb[mp_rates->irrReactions()[i]] = 0.0;
}


//==============================================================================

void Kinetics::forwardRatesOfProgress(double* const p_ropf)
{
    // Compute species concentrations (mol/m^3)
    ArrayXd conc =
        (m_thermo.numberDensity() / NA) *
        Map<const ArrayXd>(m_thermo.X(), m_thermo.nSpecies());

    forwardRatesOfProgress(conc.data(), p_ropf);
}

//==============================================================================

void Kinetics::forwardRatesOfProgress(
    const double* const p_conc, double* const p_ropf)
{
    forwardRateCoefficients(p_ropf);
    m_reactants.multReactions(p_conc, p_ropf);
    m_thirdbodies.multiplyThirdbodies(p_conc, p_ropf);
}

//==============================================================================

void Kinetics::backwardRatesOfProgress(double* const p_ropb)
{
    // Compute species concentrations (mol/m^3)
    ArrayXd conc =
        (m_thermo.numberDensity() / NA) *
        Map<const ArrayXd>(m_thermo.X(), m_thermo.nSpecies());

    backwardRatesOfProgress(conc.data(), p_ropb);
}

//==============================================================================

void Kinetics::backwardRatesOfProgress(
    const double* const p_conc, double* const p_ropb)
{
    backwardRateCoefficients(p_ropb);
    m_rev_prods.multReactions(p_conc, p_ropb);
    m_thirdbodies.multiplyThirdbodies(p_conc, p_ropb);
}

/*
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
    // Compute species concentrations (mol/m^3)
    ArrayXd conc =
        (m_thermo.numberDensity() / NA) *
        Map<const ArrayXd>(m_thermo.X(), m_thermo.nSpecies());
    
    netRatesOfProgress(conc.data(), p_rop);
}

//==============================================================================

void Kinetics::netRatesOfProgress(
    const double* const p_conc, double* const p_rop)
{
    forwardRatesOfProgress(p_conc, mp_ropf);
    backwardRatesOfProgress(p_conc, mp_ropb);

    Map<ArrayXd>(p_rop, nReactions()) = 
        Map<ArrayXd>(mp_ropf, nReactions()) - Map<ArrayXd>(mp_ropb, nReactions());
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
        std::fill(p_wdot, p_wdot + m_thermo.nSpecies(), 0);
        return;
    }

    // Compute species concentrations (mol/m^3)
    Map<ArrayXd>(p_wdot, m_thermo.nSpecies()) =
        (m_thermo.numberDensity() / NA) *
        Map<const ArrayXd>(m_thermo.X(), m_thermo.nSpecies());

    netRatesOfProgress(p_wdot, mp_rop);
    
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

void Kinetics::R(double* const p_r)
{
    // Sum up contributions from two processes:
    // diss/rec and exchange reactions.
    std::fill(p_r, p_r+m_thermo.nMolecules(), 0.0);
    R_22(mp_r_22);
    R_23(mp_r_23);

    for (int i=0; i<m_thermo.nMolecules(); ++i) 
        std::cout <<  " mp_r_22["<<i<<"] = " << mp_r_22[i] << std::endl;

    for (int i=0; i<m_thermo.nMolecules(); ++i) 
        std::cout <<  " mp_r_23["<<i<<"] = " << mp_r_23[i] << std::endl;

    // FIXME: dimensions. Should be kg / m^3-s
    // but the output from R_22 and R_23 is not ...
    for (int i=0; i<m_thermo.nMolecules(); ++i)
        p_r[i] = mp_r_22[i] + mp_r_23[i];

    for (int i=0; i<m_thermo.nMolecules(); ++i) 
        std::cout <<  " p_r["<<i<<"] = " << p_r[i] << std::endl;

    // Multiply by species molecular weights
    //for (int i = 0; i < m_thermo.nSpecies(); ++i)
    //    p_r[i] *= m_thermo.speciesMw(i);
}

//==============================================================================

// Compute the production term R_22 due to exchange reactions
void Kinetics::R_22(double* const p_r_22)
{
    double T = m_thermo.T();
    double KBT = KB * T;
    double Tv[3];
    for (int i=0; i<3; ++i)
        Tv[i] = m_thermo.Tvs(i); // NO, N2, O2

    const double * p_Y = m_thermo.Y();
    const double * p_X = m_thermo.X();
    const double   n   = m_thermo.numberDensity();

    // Partial number densities
    double nd[m_thermo.nSpecies()] = {};
    std::cout << " Number densities ... " << std::endl;
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
        nd[i] = p_X[i] * n;
	std::cout << i << " " << nd[i] << std::endl;
    }

    // Number of vibrational levels
    //int lN2 = 48; // Anharmonic
    //int lO2 = 36;
    //int lNO = 38;
    int lN2 = 33; // Harmonic
    int lO2 = 26;
    int lNO = 27;

    // Vibrational energy arrays
    static double ve[87];
    static double veN2[33];
    static double veO2[26];
    static double veNO[28];

    veNO[0] = 0.00000000e+00;
    veNO[1] = 3.78297694e-20;
    veNO[2] = 7.56595389e-20;
    veNO[3] = 1.13489308e-19;
    veNO[4] = 1.51319078e-19;
    veNO[5] = 1.89148847e-19;
    veNO[6] = 2.26978617e-19;
    veNO[7] = 2.64808386e-19;
    veNO[8] = 3.02638156e-19;
    veNO[9] = 3.40467925e-19;
    veNO[10] = 3.78297694e-19;
    veNO[11] = 4.16127464e-19;
    veNO[12] = 4.53957233e-19;
    veNO[13] = 4.91787003e-19;
    veNO[14] = 5.29616772e-19;
    veNO[15] = 5.67446542e-19;
    veNO[16] = 6.05276311e-19;
    veNO[17] = 6.43106081e-19;
    veNO[18] = 6.80935850e-19;
    veNO[19] = 7.18765620e-19;
    veNO[20] = 7.56595389e-19;
    veNO[21] = 7.94425158e-19;
    veNO[22] = 8.32254928e-19;
    veNO[23] = 8.70084697e-19;
    veNO[24] = 9.07914467e-19;
    veNO[25] = 9.45744236e-19;
    veNO[26] = 9.83574006e-19;
    veNO[27] = 1.02140378e-18;

    veN2[0] = 0.00000000e+00;
    veN2[1] = 4.68454043e-20;
    veN2[2] = 9.36908086e-20;
    veN2[3] = 1.40536213e-19;
    veN2[4] = 1.87381617e-19;
    veN2[5] = 2.34227021e-19;
    veN2[6] = 2.81072426e-19;
    veN2[7] = 3.27917830e-19;
    veN2[8] = 3.74763234e-19;
    veN2[9] = 4.21608639e-19;
    veN2[10] = 4.68454043e-19; 
    veN2[11] = 5.15299447e-19;
    veN2[12] = 5.62144851e-19;
    veN2[13] = 6.08990256e-19;
    veN2[14] = 6.55835660e-19;
    veN2[15] = 7.02681064e-19;
    veN2[16] = 7.49526469e-19;
    veN2[17] = 7.96371873e-19;
    veN2[18] = 8.43217277e-19;
    veN2[19] = 8.90062681e-19;
    veN2[20] = 9.36908086e-19;
    veN2[21] = 9.83753490e-19;
    veN2[22] = 1.03059889e-18;
    veN2[23] = 1.07744430e-18;
    veN2[24] = 1.12428970e-18;
    veN2[25] = 1.17113511e-18;
    veN2[26] = 1.21798051e-18;
    veN2[27] = 1.26482592e-18;
    veN2[28] = 1.31167132e-18;
    veN2[29] = 1.35851672e-18;
    veN2[30] = 1.40536213e-18;
    veN2[31] = 1.45220753e-18;
    veN2[32] = 1.49905294e-18;

    veO2[0] = 0.00000000e+00;
    veO2[1] = 3.13821409e-20; 
    veO2[2] = 6.27642817e-20;
    veO2[3] = 9.41464226e-20;
    veO2[4] = 1.25528563e-19;
    veO2[5] = 1.56910704e-19;
    veO2[6] = 1.88292845e-19;
    veO2[7] = 2.19674986e-19;
    veO2[8] = 2.51057127e-19;
    veO2[9] = 2.82439268e-19;
    veO2[10] = 3.13821409e-19; 
    veO2[11] = 3.45203549e-19;
    veO2[12] = 3.76585690e-19;
    veO2[13] = 4.07967831e-19;
    veO2[14] = 4.39349972e-19;
    veO2[15] = 4.70732113e-19;
    veO2[16] = 5.02114254e-19;
    veO2[17] = 5.33496395e-19;
    veO2[18] = 5.64878535e-19;
    veO2[19] = 5.96260676e-19;
    veO2[20] = 6.27642817e-19;
    veO2[21] = 6.59024958e-19;
    veO2[22] = 6.90407099e-19;
    veO2[23] = 7.21789240e-19;
    veO2[24] = 7.53171381e-19;
    veO2[25] = 7.84553521e-19;
/*
    // Number of vibrational levels
    const int lN2 = 48;
    const int lO2 = 36;
    const int lNO = 38;

    // Vibrational energy arrays in J
    static double veN2[48];
    static double veO2[36];
    static double veNO[38];

    veN2[0]  = 2.33547e-20;
    veN2[1]  = 6.96372e-20;
    veN2[2]  = 1.1535e-19 ;
    veN2[3]  = 1.60493e-19;
    veN2[4]  = 2.05065e-19;
    veN2[5]  = 2.49065e-19;
    veN2[6]  = 2.92494e-19;
    veN2[7]  = 3.35349e-19;
    veN2[8]  = 3.77629e-19;
    veN2[9]  = 4.19334e-19;
    veN2[10] = 4.60463e-19;
    veN2[11] = 5.01013e-19;
    veN2[12] = 5.40983e-19;
    veN2[13] = 5.80372e-19;
    veN2[14] = 6.19178e-19;
    veN2[15] = 6.57399e-19;
    veN2[16] = 6.95033e-19;
    veN2[17] = 7.32077e-19;
    veN2[18] = 7.68531e-19;
    veN2[19] = 8.0439e-19 ;
    veN2[20] = 8.39654e-19;
    veN2[21] = 8.74319e-19;
    veN2[22] = 9.08383e-19;
    veN2[23] = 9.41842e-19;
    veN2[24] = 9.74695e-19;
    veN2[25] = 1.00694e-18;
    veN2[26] = 1.03857e-18;
    veN2[27] = 1.06958e-18;
    veN2[28] = 1.09997e-18;
    veN2[29] = 1.12974e-18;
    veN2[30] = 1.15889e-18;
    veN2[31] = 1.1874e-18 ;
    veN2[32] = 1.21528e-18;
    veN2[33] = 1.24252e-18;
    veN2[34] = 1.26911e-18;
    veN2[35] = 1.29507e-18;
    veN2[36] = 1.32037e-18;
    veN2[37] = 1.34501e-18;
    veN2[38] = 1.369e-18  ;
    veN2[39] = 1.39232e-18;
    veN2[40] = 1.41497e-18;
    veN2[41] = 1.43695e-18;
    veN2[42] = 1.45825e-18;
    veN2[43] = 1.47887e-18;
    veN2[44] = 1.49879e-18;
    veN2[45] = 1.51803e-18;
    veN2[46] = 1.53656e-18;
    veN2[47] = 1.55439e-18;

    veO2[0]  = 1.56354e-20;
    veO2[1]  = 4.6552e-20 ;
    veO2[2]  = 7.70004e-20;
    veO2[3]  = 1.06985e-19;
    veO2[4]  = 1.3651e-19 ;
    veO2[5]  = 1.65578e-19;
    veO2[6]  = 1.94192e-19;
    veO2[7]  = 2.22354e-19;
    veO2[8]  = 2.50065e-19;
    veO2[9]  = 2.77327e-19;
    veO2[10] = 3.04138e-19;
    veO2[11] = 3.305e-19  ;
    veO2[12] = 3.56411e-19;
    veO2[13] = 3.81869e-19;
    veO2[14] = 4.06872e-19;
    veO2[15] = 4.31417e-19;
    veO2[16] = 4.55501e-19;
    veO2[17] = 4.7912e-19 ;
    veO2[18] = 5.02269e-19;
    veO2[19] = 5.24943e-19;
    veO2[20] = 5.47135e-19;
    veO2[21] = 5.6884e-19 ;
    veO2[22] = 5.90051e-19;
    veO2[23] = 6.10759e-19;
    veO2[24] = 6.30957e-19;
    veO2[25] = 6.50635e-19;
    veO2[26] = 6.69784e-19;
    veO2[27] = 6.88393e-19;
    veO2[28] = 7.06453e-19;
    veO2[29] = 7.23952e-19;
    veO2[30] = 7.40877e-19;
    veO2[31] = 7.57217e-19;
    veO2[32] = 7.72958e-19;
    veO2[33] = 7.88086e-19;
    veO2[34] = 8.02588e-19;
    veO2[35] = 8.16447e-19;

    veNO[0]  = 1.88431e-20;
    veNO[1]  = 5.61105e-20;
    veNO[2]  = 9.28207e-20;
    veNO[3]  = 1.28975e-19;
    veNO[4]  = 1.64575e-19;
    veNO[5]  = 1.99621e-19;
    veNO[6]  = 2.34116e-19;
    veNO[7]  = 2.68059e-19;
    veNO[8]  = 3.01454e-19;
    veNO[9]  = 3.343e-19  ;
    veNO[10] = 3.666e-19  ;
    veNO[11] = 3.98354e-19;
    veNO[12] = 4.29564e-19;
    veNO[13] = 4.60232e-19;
    veNO[14] = 4.90357e-19;
    veNO[15] = 5.19943e-19;
    veNO[16] = 5.4899e-19 ;
    veNO[17] = 5.77499e-19;
    veNO[18] = 6.05472e-19;
    veNO[19] = 6.3291e-19 ;
    veNO[20] = 6.59815e-19;
    veNO[21] = 6.86187e-19;
    veNO[22] = 7.12028e-19;
    veNO[23] = 7.3734e-19 ;
    veNO[24] = 7.62123e-19;
    veNO[25] = 7.86379e-19;
    veNO[26] = 8.10109e-19;
    veNO[27] = 8.33315e-19;
    veNO[28] = 8.55998e-19;
    veNO[29] = 8.78159e-19;
    veNO[30] = 8.99799e-19;
    veNO[31] = 9.2092e-19 ;
    veNO[32] = 9.41523e-19;
    veNO[33] = 9.6161e-19 ;
    veNO[34] = 9.81182e-19;
    veNO[35] = 1.00024e-18;
    veNO[36] = 1.01878e-18;
    veNO[37] = 1.03682e-18;
*/
    // Dissociation energy (ground level) in J
    const double de[m_thermo.nMolecules()] = {1.0409017238088574e-18,
        1.5636156480913654e-18, 8.196091362099268e-19};

    // Non-equilibrium partition functions (anharmonic)
    double z_vibr_T_Tv_NO = 0.;
    for (int i=0; i<lNO; ++i) {

	// anharmonic (treanor-marrone)
        //z_vibr_T_Tv_NO += exp(-(veNO[i]-i*veNO[0])/KBT-
	//		                i*veNO[0]/(KB*Tv[0]));
	// harmonic (boltzmann)
        z_vibr_T_Tv_NO += exp(-veNO[i]/(KB*Tv[0]));
    }

    double z_vibr_T_Tv_N2 = 0.;
    for (int i=0; i<lN2; ++i) {
        //z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/KBT-
	//		                i*veN2[0]/(KB*Tv[1]));
        z_vibr_T_Tv_N2 += exp(-veN2[i]/(KB*Tv[1]));
    }

    double z_vibr_T_Tv_O2 = 0.;
    for (int i=0; i<lO2; ++i) {
        //z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/KBT-
	//		                i*veO2[0]/(KB*Tv[2]));
        z_vibr_T_Tv_O2 += exp(-veO2[i]/(KB*Tv[2]));
    }

    std::cout << "Non-equilibrium partition functions ... " << std::endl;
    std::cout << z_vibr_T_Tv_NO << " " 
	      << z_vibr_T_Tv_N2 << " " 	
	      << z_vibr_T_Tv_O2 << std::endl;

    // Compute the vibrational populations according to 
    // Treanor-Marrone or Boltzmann distributions
    double nNO[lNO];
    double nN2[lN2];
    double nO2[lO2];

    std::cout << " Vibrational populations ... " << std::endl;
    for (int i=0; i<lNO; i++) {
        //nNO[i] = nd[2] * exp(-(veNO[i]-i*veNO[0])/(KB*T)+i*veNO[0]/(KB*Tv[0]))
	//    / z_vibr_T_Tv_NO;

	// Boltzmann distribution
        nNO[i] = nd[2] * exp(-veNO[i]/(KB*Tv[0])) / z_vibr_T_Tv_NO;
	std::cout << i << " " << nNO[i] << std::endl;
    }

    for (int i=0; i<lN2; i++) {
        //nN2[i] = nd[3] * exp(-(veN2[i]-i*veN2[0])/(KB*T)+i*veN2[0]/(KB*Tv[1]))
	//    / z_vibr_T_Tv_N2;
        nN2[i] = nd[3] * exp(-veN2[i]/(KB*Tv[1])) / z_vibr_T_Tv_N2;
	std::cout << i << " " << nN2[i] << std::endl;
    }

    for (int i=0; i<lO2; i++) {
        //nO2[i] = nd[4] * exp(-(veO2[i]-i*veO2[0])/(KB*T)+i*veO2[0]/(KB*Tv[2]))
	//    / z_vibr_T_Tv_O2;
        nO2[i] = nd[4] * exp(-veO2[i]/(KB*Tv[2])) / z_vibr_T_Tv_O2;
	std::cout << i << " " << nO2[i] << std::endl;
    }

    // Mass and diameter
    double m[5], d[5]; // kg, m
    m[0] = 2.3258672171567998e-26; // N
    m[1] = 2.6567628316576e-26;    // O
    m[2] = 4.9826300488143997e-26; // NO
    m[3] = 4.6517344343135997e-26; // N2
    m[4] = 5.3135256633152e-26;	   // O2
    d[0] = 3.298e-10;
    d[1] = 2.75e-10;
    d[2] = 3.4061e-10;
    d[3] = 3.40385035355259e-10;
    d[4] = 3.5155e-10;

    // Computation of STS rate coefficients for the exchange reactions
    // with the model of Savelev. We only consider 2 reactions:
    // N2(i) + O <=> NO(k) + N
    // O2(i) + N <=> NO(k) + O
    double k_N2_O[lN2][lNO];
    double k_O2_N[lO2][lNO];

    // constants for N2,O2
    double A[2] = {8e-17, 4e-16}; // m^3/sec
    double b[2] = {0, -0.39};
    double U = 3.*T; // A/(6.*KB);
    double KBU = KB * U;

    // N2 vibr. energy, J
    double en2[lN2]    = {};
    double en2_eV[lN2] = {};
    for (int i=0; i<lN2; ++i) {
        en2[i] = veN2[i] + veN2[0];    // J
        en2_eV[i] = en2[i] * 6.242e18; // eV
	std::cout << " en2_eV: " << i << " " << en2_eV[i] << std::endl;
    }

    // O2 vibr. energy, J
    double eo2[lO2]    = {};
    double eo2_eV[lO2] = {};
    for (int i=0; i<lO2; ++i) {
        eo2[i] = veO2[i] + veO2[0];
        eo2_eV[i] = eo2[i] * 6.242e18; 
	std::cout << " eo2_eV: " << i << " " << eo2_eV[i] << std::endl;
    }

    double eno[lNO]      = {};
    double eno_eV[lNO]   = {};
    double Ear_n2[lNO]   = {};
    double Ear_n2_J[lNO] = {};
    double Ear_o2[lNO]   = {};
    double Ear_o2_J[lNO] = {};
    for (int i=0; i<lNO; ++i) {
        eno[i] = veNO[i] + veNO[0];
        eno_eV[i] = eno[i] * 6.242e18; 
	std::cout << " eno_eV: " << i << " " << eno_eV[i] << std::endl;
    }

    std::cout << " Ear ... " << std::endl;
    for (int i=0; i<lNO; ++i) {

	// N2+O reaction
        Ear_n2[i] = 2.8793 + 1.02227 * eno_eV[i]; // eV
        Ear_n2_J[i] = Ear_n2[i] / 6.242e18;       // J

	// O2+N reaction
	if (eno_eV[i] < 1.3706) {
            Ear_o2[i] = 0.098;
	}
	else if ((1.3706 < eno_eV[i]) && (eno_eV[i] < 2.4121)) {
	    Ear_o2[i] = -0.6521 + 0.54736 * eno_eV[i];
	}
	else if (eno_eV[i] > 2.4121) {
            Ear_o2[i] = -1.8451 + 1.04189 * eno_eV[i]; // eV
	} else {
            std::cout << "Something wrong in Ear_o2" << std::endl;
	}
        Ear_o2_J[i] = Ear_o2[i] / 6.242e18; // J

	std::cout << i << " " << Ear_n2_J[i] << "  "
	 	  	      << Ear_o2_J[i] << std::endl;
    }

    // Vibration partition functions
    std::cout << "Vibr. partition functions" << std::endl;
    double Zv_n2;
    for (int i=0; i<lN2; ++i)
        Zv_n2 += exp(-en2[i]/KBT);

    double Zv_o2;
    for (int i=0; i<lO2; ++i) 
        Zv_o2 += exp(-eo2[i]/KBT);

    std::cout << " Zv_n2: " << Zv_n2 << " " 
	      << " Zv_o2: " << Zv_o2 << std::endl;

    // Equilibrium coefficients
    std::cout << " N2 eq. coeffs..." << std::endl;
    double k_eq_n2[lNO] = {}; // m^3/sec
    double fac = A[0] * pow(T,b[0]);
    for (int i=0; i<lNO; ++i) {
        k_eq_n2[i] = fac * (1. + 0.333333333333333 * eno_eV[i]) 
	    * exp(-Ear_n2_J[i] / KBT); 
	std::cout << i << " " << k_eq_n2[i] << std::endl;
    }

    std::cout << " O2 eq. coeffs..." << std::endl;
    double k_eq_o2[lNO] = {};
    fac = A[1] * pow(T,b[1]);
    for (int i=0; i<lNO; ++i) {
        k_eq_o2[i] = fac * (Ear_o2[i] + 0.8) * exp(-Ear_o2_J[i] / KBT);
        std::cout << i << " " << k_eq_o2[i] << std::endl;
    }

    // Energy threshold, for each k -> e_i*
    std::cout << "energy threshold" << std::endl;
    double sum1_n2[lNO] = {}; //zeros(lno,1);
    double sum2_n2[lNO] = {}; 
    double sum1_o2[lNO] = {};
    double sum2_o2[lNO] = {};

    // Identify the thresholds to compute the partial sums
    int thrsN2[lNO];
    int thrsO2[lNO];
    for (int i=0; i<lNO; ++i) {

        for (int j=0; j<lN2; ++j) {
	    if (en2[j] <= Ear_n2_J[i]) thrsN2[i] = j;
	}

        for (int j=0; j<lO2; ++j) {
	    if (eo2[j] <= Ear_o2_J[i]) thrsO2[i] = j;
	}
        std::cout << " thrsN2["<<i<<"] = " << thrsN2[i] 
		  << " thrsO2["<<i<<"] = " << thrsO2[i] << std::endl;
    }

    std::cout << " Partial sums " << std::endl;
    for (int i=0; i<lNO; ++i) {

        // N2 + O reaction
        for (int j=0; j<thrsN2[i]; ++j) {
	    sum1_n2[i] += exp(-(Ear_n2_J[i]-en2[j])/KBU);    
	}
        for (int j=thrsN2[i]; j<lN2; ++j) {
	    sum2_n2[i] += exp(-(Ear_n2_J[i]-en2[j])/KBT);    
	}

        // O2 + N reaction
        for (int j=0; j<thrsO2[i]; ++j) {
	    sum1_o2[i] += exp(-(Ear_o2_J[i]-eo2[j])/KBU);    
	}
        for (int j=thrsN2[i]; j<lO2; ++j) {
	    sum2_o2[i] += exp(-(Ear_o2_J[i]-eo2[j])/KBT);    
	}
        std::cout << i << " " << sum1_n2[i] << " " << sum2_n2[i] << " " 
        	              << sum1_o2[i] << " " << sum2_o2[i] << std::endl;
    }

    // Normalizing coefficients
    std::cout << "Normalizing coeffs" << std::endl;
    double C_n2[lNO]  = {}, C_o2[lNO] = {};
    double B1_n2[lNO] = {}, B2_n2[lNO] = {};
    double B1_o2[lNO] = {}, B2_o2[lNO] = {}; 

    for (int i=0; i<lNO; ++i) {
        C_n2[i] = Zv_n2 * 1./(sum1_n2[i] + sum2_n2[i]);
        B1_n2[i] = C_n2[i] * k_eq_n2[i] * exp(-Ear_n2_J[i]/KBU);
        B2_n2[i] = C_n2[i] * k_eq_n2[i] * exp(Ear_n2_J[i]/KBT);
        std::cout << " C_n2["<<i<<"] = " << C_n2[i] 
		  << " B1_n2["<<i<<"] = " << B1_n2[i] 
		  << " B2_n2["<<i<<"] = " << B2_n2[i] << std::endl;
    }

    for (int i=0; i<lNO; ++i) {
        C_o2[i] = Zv_o2 * 1./(sum1_o2[i] + sum2_o2[i]);
        B1_o2[i] = C_o2[i] * k_eq_o2[i] * exp(-Ear_o2_J[i]/KBU);
        B2_o2[i] = C_o2[i] * k_eq_o2[i] * exp(Ear_o2_J[i]/KBT);
        std::cout << " C_o2["<<i<<"] = " << C_o2[i] 
		  << " B1_o2["<<i<<"] = " << B1_o2[i] 
		  << " B2_o2["<<i<<"] = " << B2_o2[i] << std::endl;
    }

    // TODO: fix bug in the exponential
    std::cout << " k_N2_O " << std::endl;
    double kbTU =  (1./T + 1./U) / KB;
    for (int i=0; i<lN2; ++i) { // N2 rows
        for (int j = 0; j < lNO; ++j) { // NO columns
	    if (en2[i] < Ear_n2_J[j]) { // Ear_n2_J has size = lNO 
	        k_N2_O[i][j] = B1_n2[i] * exp(en2[i]/kbTU);
	    }
	    else if (en2[i] > Ear_n2_J[j]) {
		k_N2_O[i][j] = B2_n2[i];
	    } else {
                std::cout << " Out-of-bound error ... " << std::endl;
	    }
            std::cout << " k_N2_O["<<i<<"]["<<j<<"] = " 
		      << k_N2_O[i][j] << std::endl;
	}
    }

    // TODO: fix bug in the exponential
    std::cout << " k_O2_N " << std::endl;
    for (int i=0; i<lO2; ++i) { // O2 rows 
        for (int j = 0; j < lNO; ++j) { // NO columns
	    if (eo2[i] < Ear_o2_J[j]) { // Ear_o2_J has size = lNO
	        k_O2_N[i][j] = B1_o2[i] * exp(eo2[i]/kbTU);
	    }
	    else if (eo2[i] > Ear_o2_J[j]) {
		k_O2_N[i][j] = B2_o2[i];
	    } else {
                std::cout << " Out-of-bound error ... " << std::endl;
	    }
            std::cout << " k_O2_N["<<i<<"]["<<j<<"] = " 
		      << k_O2_N[i][j] << std::endl;
	}
    }

    /*
    ////////////////////////////////////////////////////////////////////////////
    // 1. Rusanov-Friedman model
    for (int i=0; i<lN2; ++i) {
        if (Ear_n2_J[i] > 0.51 * en2[i] ) {
          B2_n2[i] = A[0]*pow(T, b[0])*exp(-(Ear_n2_J[i] - 0.51*en2[i])/(KB*T));
        } else {
          B2_n2[i] = A[0] * pow(T, b[0]);
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    // 2. Polak model
    for (int i=0; i<lN2; ++i) {
        if (Ear_n2_J[i] > 0.52 * en2[i] ) {
          B2_n2[i] = A[0] * pow(T, b[0]) *
	      exp(-(Ear_n2_J[i] - 0.52 * en2[i]) / (0.9 * KB * T));
        } else {
          B2_n2[i] = A[0] * pow(T, b[0]);
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    // 3. Warnatz model
    for (int i=0; i<lN2; ++i) {
        if (Ear_n2_J[i] > en2[i] ) {
          B2_n2[i] = (4.17e+12 + i) * pow(T, 0) *
	      exp(-(Ear_n2_J[i] - en2[i]) / (KB * T));
        } else {
          B2_n2[i] = (4.17e+12 + i) * pow(T, 0);
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    */

    // Compute the rotational partition function
    // according to a simplified formula ... TODO: better formula
    double Be[3] = {167.20, 199.8, 143.77}; // m^-1
    int sigma[3] = {1, 2, 2};
    double theta_r[3] = {};
    double z_rot[3] = {};
    for (int i=0; i<3; ++i) {
        theta_r[i] = Be[i]*HP*C0/KB;
        z_rot[i] = T/(sigma[i]*theta_r[i]);
    }

    std::cout << " Kz_n2 ... " << std::endl;
    double Kz_n2[lN2][lNO] = {}; 
    double fac0 = pow((m[3]*m[1]/(m[2]*m[0])),1.5) * z_rot[1]/z_rot[0] 
	            * exp((de[1]-de[0])/KBT);
    for (int i=0; i<lN2; ++i) {
        for (int j=0; j<lNO; ++j) {
            Kz_n2[i][j] = fac0 * exp((-veN2[i] + veNO[j])/KBT);
	    cout << Kz_n2[i][j] << " \n"[j == lNO-1];
	}
    }

    std::cout << "Kz_o2" << std::endl;
    double Kz_o2[lO2][lNO] = {}; 
    double fac1 = pow((m[4]*m[0]/(m[2]*m[1])),1.5) * z_rot[2]/z_rot[0] 
	            * exp((de[2]-de[0])/KBT);
    for (int i=0; i<lO2; ++i) {
        for (int j=0; j<lNO; ++j) {
            Kz_o2[i][j] = fac1 * exp((-veNO[j] + veO2[i])/KBT);
	    cout << Kz_o2[i][j] << " \n"[j == lNO-1];
	}
    }

    // Compute MT forward exchange coefficients, eq. 3.93.
    // N2
    double k_O_N2[lN2][lNO];
    double fac0f, fac0b = 0.;
    for (int i=0; i<lN2; ++i) {
        for (int j=0; j<lNO; ++j) {

	    // MT forward rates
            fac0f += exp(-((veN2[i]-i*veN2[0])/KBT) -
	                            i*veN2[0]/(KB*Tv[1])) * k_N2_O[i][j];
	    // STS backward rates
	    k_O_N2[i][j] = Kz_n2[i][j] * k_N2_O[i][j];
            fac0b += exp(-((veN2[i]-i*veN2[0])/KBT) -
	                            i*veN2[0]/(KB*Tv[1])) * k_O_N2[i][j];
	}
    }
    double k_exch_N2_NO = fac0f / z_vibr_T_Tv_N2;
    double k_exch_NO_N2 = fac0b / z_vibr_T_Tv_N2;
    std::cout << " k_exch_N2_NO: " << k_exch_N2_NO << std::endl; 
    std::cout << " k_exch_NO_N2: " << k_exch_NO_N2 << std::endl; 

    // O2
    double k_N_O2[lO2][lNO];
    for (int i=0; i<lO2; ++i) {
	fac0f = fac0b = 0.;
        for (int j=0; j<lNO; ++j) {

	    // MT forward rates
            fac0f += exp(-((veO2[i]-i*veO2[0])/KBT) -
                                    i*veO2[0]/(KB*Tv[2])) * k_O2_N[i][j];
	    // STS backward rates
	    k_N_O2[i][j] = Kz_o2[i][j] * k_N_O2[i][j];
            fac0b += exp(-((veO2[i]-i*veO2[0])/KBT) -
                                    i*veO2[0]/(KB*Tv[2])) * k_N_O2[i][j];
        }
    }
    double k_exch_O2_NO = fac0f / z_vibr_T_Tv_O2;
    double k_exch_NO_O2 = fac0b / z_vibr_T_Tv_O2;
    std::cout << " k_exch_O2_NO: " << k_exch_O2_NO << std::endl; 
    std::cout << " k_exch_NO_O2: " << k_exch_NO_O2 << std::endl; 

    double R_exch_22[m_thermo.nMolecules()] = {};

    // NO
    R_exch_22[0] = (nd[3] * nd[1] * k_exch_N2_NO - nd[2] * nd[0] * k_exch_NO_N2) + 
	           (nd[4] * nd[0] * k_exch_O2_NO - nd[2] * nd[1] * k_exch_NO_O2);
    // N2
    R_exch_22[1] = nd[2] * nd[0] * k_exch_NO_N2 - nd[3] * nd[1] * k_exch_N2_NO;
    // O2
    R_exch_22[2] = nd[2] * nd[1] * k_exch_NO_O2 - nd[4] * nd[0] * k_exch_O2_NO;

    // Compute species concentrations (mol/m^3)
    Map<ArrayXd>(p_r_22, m_thermo.nMolecules()) = 0.;
    
    // FIXME: dimensions ...
    std::cout << " R_exch_22[i] " << std::endl;
    for (int i=0; i<m_thermo.nMolecules(); ++i) {
        p_r_22[i] = R_exch_22[i];
	std::cout << i << " " << p_r_22[i] << std::endl;
    }
}

//==============================================================================

  // The influence of vibrational state-resolvent transport coefficients on the 
  // wave propagation in diatomic gases.
  // Compute rate coefficients of the VT transition using 
  // Schwartz-Slavsky-Herzfeld (SSH) theory for the harmonic oscillator model
  double Kinetics::k_VT_SSH(double T, int i, double coll_mass, double diameter, 
      double omega_e, double epsilon, double r_e) 
  {
    // eq. 80b
    return p_Z_coll(T, 1., coll_mass, diameter) * i * 
        P_SSH_VT_10(T, coll_mass, omega_e, epsilon, diameter, r_e);
  }

//==============================================================================

  double Kinetics::k_VT_SSH(double T, int i, double coll_mass, double diameter,
      double omega_e, double epsilon, double r_e, double Delta_E_vibr,
          double vibr_energy_1) 
  {
      double alpha = 17.5 / diameter;
      double gamma0 = PI * sqrt(coll_mass / (2 * KB * T)) / (alpha * HPBAR);
      double gammai = gamma0 * (vibr_energy_1 - 2 * i * Delta_E_vibr);
      gamma0 *= vibr_energy_1;
      double delta_VT = Delta_E_vibr / vibr_energy_1;

      if (gammai >= 20) {
        delta_VT *= 4 * pow(gamma0, 2./3.);
      } else {
        delta_VT *= (4./3) * gamma0;
      }

      return k_VT_SSH(T, i, coll_mass, diameter, omega_e, epsilon, r_e) * 
          exp(i * delta_VT) * exp(-i * Delta_E_vibr / (KB * T));
  }

//==============================================================================
 
  // Compute the average probability of the VT transition M(1) + M -> M(0) + M
  double Kinetics::P_SSH_VT_10(double T, double coll_mass, double omega_e, 
      double epsilon, double diameter, double r_e) 
  {
      double kT = KB * T;
      // eq. 86a
      double alpha = 17.5 / diameter;
      double omega = 2 * PI * C0 * omega_e;
      // eq. 86b
      double chi = pow(PI * PI * coll_mass * omega * omega / 
          (2 * alpha * alpha * kT), 0.3333333); // TODO 4 or 2?
      // eq. 86c
      double r = diameter * 
          pow(0.5 * sqrt(1 + chi * kT / epsilon) + 0.5, -0.1666666);
      // eq. 85
      double Z0 = alpha * r_e * alpha * r_e * 
          exp( -(3 * alpha * r_e * r_e) / (8 * r) );
      // eq. 83
      return 1.294 * pow(r / diameter, 2) / (Z0 * alpha * alpha * HPBAR * 
          (1 + 1.1 * epsilon / kT)) * 4 * PI * PI * coll_mass * omega * 
	      sqrt(4 * PI * chi / 3) * 
	          exp(-3*chi + (HPBAR * omega / (2 * kT)) + epsilon/kT);
  }

//==============================================================================

  double Kinetics::p_Z_coll(double T, double n, double coll_mass, 
      double diameter) 
  {
      // collision frequency, eq. 82a
      return 4. * PI * n * diameter * diameter
                * sqrt(KB * T / (2. * PI * coll_mass));
  }

//==============================================================================

// Compute the production term R_23 due to dissociation/recombination
void Kinetics::R_23(double* const p_r_23)
{
    //const double * p_Y = m_thermo.Y();
    const double * p_X = m_thermo.X();
    const double   n   = m_thermo.numberDensity();
    std::cout << "number density: " << n << std::endl;
    double T = m_thermo.T();
    std::cout << T << std::endl;
    double Tv[3];
    for (int i=0; i<3; ++i) {
        Tv[i] = m_thermo.Tvs(i); // NO, N2, O2
        std::cout << Tv[i] << std::endl;
    }
    double KBT = KB * T;

    // Number of vibrational levels
    //int lN2 = 48; // Anharmonic
    //int lO2 = 36;
    //int lNO = 38;
    int lN2 = 33; // Harmonic
    int lO2 = 26;
    int lNO = 27;

    // Vibrational energy arrays
    //static double ve[122];
    //static double veN2[48];
    //static double veO2[36];
    //static double veNO[38];
    static double ve[87];
    static double veN2[33];
    static double veO2[26];
    static double veNO[28];

//    ve[1] = 6.96372e-20;
//    ve[2] = 1.1535e-19 ;
//    ve[3] = 1.60493e-19;
//    ve[4] = 2.05065e-19;
//    ve[5] = 2.49065e-19;
//    ve[6] = 2.92494e-19;
//    ve[7] = 3.35349e-19;
//    ve[8] = 3.77629e-19;
//    ve[9] = 4.19334e-19;
//    ve[10] = 4.60463e-19;
//    ve[11] = 5.01013e-19;
//    ve[12] = 5.40983e-19;
//    ve[13] = 5.80372e-19;
//    ve[14] = 6.19178e-19;
//    ve[15] = 6.57399e-19;
//    ve[16] = 6.95033e-19;
//    ve[17] = 7.32077e-19;
//    ve[18] = 7.68531e-19;
//    ve[19] = 8.0439e-19;
//    ve[20] = 8.39654e-19;
//    ve[21] = 8.74319e-19;
//    ve[22] = 9.08383e-19;
//    ve[23] = 9.41842e-19;
//    ve[24] = 9.74695e-19;
//    ve[25] = 1.00694e-18;
//    ve[26] = 1.03857e-18;
//    ve[27] = 1.06958e-18;
//    ve[28] = 1.09997e-18;
//    ve[29] = 1.12974e-18;
//    ve[30] = 1.15889e-18;
//    ve[31] = 1.1874e-18;
//    ve[32] = 1.21528e-18;
//    ve[33] = 1.24252e-18;
//    ve[34] = 1.26911e-18;
//    ve[35] = 1.29507e-18;
//    ve[36] = 1.32037e-18;
//    ve[37] = 1.34501e-18;
//    ve[38] = 1.369e-18  ;
//    ve[39] = 1.39232e-18;
//    ve[40] = 1.41497e-18;
//    ve[41] = 1.43695e-18;
//    ve[42] = 1.45825e-18;
//    ve[43] = 1.47887e-18;
//    ve[44] = 1.49879e-18;
//    ve[45] = 1.51803e-18;
//    ve[46] = 1.53656e-18;
//    ve[47] = 1.55439e-18;
//    ve[48] = 1.56354e-20;
//    ve[49] = 4.6552e-20;
//    ve[50] = 7.70004e-20;
//    ve[51] = 1.06985e-19;
//    ve[52] = 1.3651e-19;
//    ve[53] = 1.65578e-19;
//    ve[54] = 1.94192e-19;
//    ve[55] = 2.22354e-19;
//    ve[56] = 2.50065e-19;
//    ve[57] = 2.77327e-19;
//    ve[58] = 3.04138e-19;
//    ve[59] = 3.305e-19  ;
//    ve[60] = 3.56411e-19;
//    ve[61] = 3.81869e-19;
//    ve[62] = 4.06872e-19;
//    ve[63] = 4.31417e-19;
//    ve[64] = 4.55501e-19;
//    ve[65] = 4.7912e-19;
//    ve[66] = 5.02269e-19;
//    ve[67] = 5.24943e-19;
//    ve[68] = 5.47135e-19;
//    ve[69] = 5.6884e-19;
//    ve[70] = 5.90051e-19;
//    ve[71] = 6.10759e-19;
//    ve[72] = 6.30957e-19;
//    ve[73] = 6.50635e-19;
//    ve[74] = 6.69784e-19;
//    ve[75] = 6.88393e-19;
//    ve[76] = 7.06453e-19;
//    ve[77] = 7.23952e-19;
//    ve[78] = 7.40877e-19;
//    ve[79] = 7.57217e-19;
//    ve[80] = 7.72958e-19;
//    ve[81] = 7.88086e-19;
//    ve[82] = 8.02588e-19;
//    ve[83] = 8.16447e-19;
//    ve[84] = 1.88431e-20;
//    ve[85] = 5.61105e-20;
//    ve[86] = 9.28207e-20;
//    ve[87] = 1.28975e-19;
//    ve[88] = 1.64575e-19;
//    ve[89] = 1.99621e-19;
//    ve[90] = 2.34116e-19;
//    ve[91] = 2.68059e-19;
//    ve[92] = 3.01454e-19;
//    ve[93] = 3.343e-19  ;
//    ve[94] = 3.666e-19  ;
//    ve[95] = 3.98354e-19;
//    ve[96] = 4.29564e-19;
//    ve[97] = 4.60232e-19;
//    ve[98] = 4.90357e-19;
//    ve[99] = 5.19943e-19;
//    ve[100] = 5.4899e-19;
//    ve[101] = 5.77499e-19;
//    ve[102] = 6.05472e-19;
//    ve[103] = 6.3291e-19 ;
//    ve[104] = 6.59815e-19;
//    ve[105] = 6.86187e-19;
//    ve[106] = 7.12028e-19;
//    ve[107] = 7.3734e-19 ;
//    ve[108] = 7.62123e-19;
//    ve[109] = 7.86379e-19;
//    ve[110] = 8.10109e-19;
//    ve[111] = 8.33315e-19;
//    ve[112] = 8.55998e-19;
//    ve[113] = 8.78159e-19;
//    ve[114] = 8.99799e-19;
//    ve[115] = 9.2092e-19 ;
//    ve[116] = 9.41523e-19;
//    ve[117] = 9.6161e-19 ;
//    ve[118] = 9.81182e-19;
//    ve[119] = 1.00024e-18;
//    ve[120] = 1.01878e-18;
//    ve[121] = 1.03682e-18;
//
//    veN2[0]  = 2.33547e-20;
//    veN2[1]  = 6.96372e-20;
//    veN2[2]  = 1.1535e-19 ;
//    veN2[3]  = 1.60493e-19;
//    veN2[4]  = 2.05065e-19;
//    veN2[5]  = 2.49065e-19;
//    veN2[6]  = 2.92494e-19;
//    veN2[7]  = 3.35349e-19;
//    veN2[8]  = 3.77629e-19;
//    veN2[9]  = 4.19334e-19;
//    veN2[10] = 4.60463e-19;
//    veN2[11] = 5.01013e-19;
//    veN2[12] = 5.40983e-19;
//    veN2[13] = 5.80372e-19;
//    veN2[14] = 6.19178e-19;
//    veN2[15] = 6.57399e-19;
//    veN2[16] = 6.95033e-19;
//    veN2[17] = 7.32077e-19;
//    veN2[18] = 7.68531e-19;
//    veN2[19] = 8.0439e-19 ;
//    veN2[20] = 8.39654e-19;
//    veN2[21] = 8.74319e-19;
//    veN2[22] = 9.08383e-19;
//    veN2[23] = 9.41842e-19;
//    veN2[24] = 9.74695e-19;
//    veN2[25] = 1.00694e-18;
//    veN2[26] = 1.03857e-18;
//    veN2[27] = 1.06958e-18;
//    veN2[28] = 1.09997e-18;
//    veN2[29] = 1.12974e-18;
//    veN2[30] = 1.15889e-18;
//    veN2[31] = 1.1874e-18 ;
//    veN2[32] = 1.21528e-18;
//    veN2[33] = 1.24252e-18;
//    veN2[34] = 1.26911e-18;
//    veN2[35] = 1.29507e-18;
//    veN2[36] = 1.32037e-18;
//    veN2[37] = 1.34501e-18;
//    veN2[38] = 1.369e-18  ;
//    veN2[39] = 1.39232e-18;
//    veN2[40] = 1.41497e-18;
//    veN2[41] = 1.43695e-18;
//    veN2[42] = 1.45825e-18;
//    veN2[43] = 1.47887e-18;
//    veN2[44] = 1.49879e-18;
//    veN2[45] = 1.51803e-18;
//    veN2[46] = 1.53656e-18;
//    veN2[47] = 1.55439e-18;
//
//    veO2[0]  = 1.56354e-20;
//    veO2[1]  = 4.6552e-20 ;
//    veO2[2]  = 7.70004e-20;
//    veO2[3]  = 1.06985e-19;
//    veO2[4]  = 1.3651e-19 ;
//    veO2[5]  = 1.65578e-19;
//    veO2[6]  = 1.94192e-19;
//    veO2[7]  = 2.22354e-19;
//    veO2[8]  = 2.50065e-19;
//    veO2[9]  = 2.77327e-19;
//    veO2[10] = 3.04138e-19;
//    veO2[11] = 3.305e-19  ;
//    veO2[12] = 3.56411e-19;
//    veO2[13] = 3.81869e-19;
//    veO2[14] = 4.06872e-19;
//    veO2[15] = 4.31417e-19;
//    veO2[16] = 4.55501e-19;
//    veO2[17] = 4.7912e-19 ;
//    veO2[18] = 5.02269e-19;
//    veO2[19] = 5.24943e-19;
//    veO2[20] = 5.47135e-19;
//    veO2[21] = 5.6884e-19 ;
//    veO2[22] = 5.90051e-19;
//    veO2[23] = 6.10759e-19;
//    veO2[24] = 6.30957e-19;
//    veO2[25] = 6.50635e-19;
//    veO2[26] = 6.69784e-19;
//    veO2[27] = 6.88393e-19;
//    veO2[28] = 7.06453e-19;
//    veO2[29] = 7.23952e-19;
//    veO2[30] = 7.40877e-19;
//    veO2[31] = 7.57217e-19;
//    veO2[32] = 7.72958e-19;
//    veO2[33] = 7.88086e-19;
//    veO2[34] = 8.02588e-19;
//    veO2[35] = 8.16447e-19;
//
//    veNO[0]  = 1.88431e-20;
//    veNO[1]  = 5.61105e-20;
//    veNO[2]  = 9.28207e-20;
//    veNO[3]  = 1.28975e-19;
//    veNO[4]  = 1.64575e-19;
//    veNO[5]  = 1.99621e-19;
//    veNO[6]  = 2.34116e-19;
//    veNO[7]  = 2.68059e-19;
//    veNO[8]  = 3.01454e-19;
//    veNO[9]  = 3.343e-19  ;
//    veNO[10] = 3.666e-19  ;
//    veNO[11] = 3.98354e-19;
//    veNO[12] = 4.29564e-19;
//    veNO[13] = 4.60232e-19;
//    veNO[14] = 4.90357e-19;
//    veNO[15] = 5.19943e-19;
//    veNO[16] = 5.4899e-19 ;
//    veNO[17] = 5.77499e-19;
//    veNO[18] = 6.05472e-19;
//    veNO[19] = 6.3291e-19 ;
//    veNO[20] = 6.59815e-19;
//    veNO[21] = 6.86187e-19;
//    veNO[22] = 7.12028e-19;
//    veNO[23] = 7.3734e-19 ;
//    veNO[24] = 7.62123e-19;
//    veNO[25] = 7.86379e-19;
//    veNO[26] = 8.10109e-19;
//    veNO[27] = 8.33315e-19;
//    veNO[28] = 8.55998e-19;
//    veNO[29] = 8.78159e-19;
//    veNO[30] = 8.99799e-19;
//    veNO[31] = 9.2092e-19 ;
//    veNO[32] = 9.41523e-19;
//    veNO[33] = 9.6161e-19 ;
//    veNO[34] = 9.81182e-19;
//    veNO[35] = 1.00024e-18;
//    veNO[36] = 1.01878e-18;
//    veNO[37] = 1.03682e-18;

    veNO[0] = 0.00000000e+00;
    veNO[1] = 3.78297694e-20;
    veNO[2] = 7.56595389e-20;
    veNO[3] = 1.13489308e-19;
    veNO[4] = 1.51319078e-19;
    veNO[5] = 1.89148847e-19;
    veNO[6] = 2.26978617e-19;
    veNO[7] = 2.64808386e-19;
    veNO[8] = 3.02638156e-19;
    veNO[9] = 3.40467925e-19;
    veNO[10] = 3.78297694e-19;
    veNO[11] = 4.16127464e-19;
    veNO[12] = 4.53957233e-19;
    veNO[13] = 4.91787003e-19;
    veNO[14] = 5.29616772e-19;
    veNO[15] = 5.67446542e-19;
    veNO[16] = 6.05276311e-19;
    veNO[17] = 6.43106081e-19;
    veNO[18] = 6.80935850e-19;
    veNO[19] = 7.18765620e-19;
    veNO[20] = 7.56595389e-19;
    veNO[21] = 7.94425158e-19;
    veNO[22] = 8.32254928e-19;
    veNO[23] = 8.70084697e-19;
    veNO[24] = 9.07914467e-19;
    veNO[25] = 9.45744236e-19;
    veNO[26] = 9.83574006e-19;
    veNO[27] = 1.02140378e-18;

    veN2[0] = 0.00000000e+00;
    veN2[1] = 4.68454043e-20;
    veN2[2] = 9.36908086e-20;
    veN2[3] = 1.40536213e-19;
    veN2[4] = 1.87381617e-19;
    veN2[5] = 2.34227021e-19;
    veN2[6] = 2.81072426e-19;
    veN2[7] = 3.27917830e-19;
    veN2[8] = 3.74763234e-19;
    veN2[9] = 4.21608639e-19;
    veN2[10] = 4.68454043e-19; 
    veN2[11] = 5.15299447e-19;
    veN2[12] = 5.62144851e-19;
    veN2[13] = 6.08990256e-19;
    veN2[14] = 6.55835660e-19;
    veN2[15] = 7.02681064e-19;
    veN2[16] = 7.49526469e-19;
    veN2[17] = 7.96371873e-19;
    veN2[18] = 8.43217277e-19;
    veN2[19] = 8.90062681e-19;
    veN2[20] = 9.36908086e-19;
    veN2[21] = 9.83753490e-19;
    veN2[22] = 1.03059889e-18;
    veN2[23] = 1.07744430e-18;
    veN2[24] = 1.12428970e-18;
    veN2[25] = 1.17113511e-18;
    veN2[26] = 1.21798051e-18;
    veN2[27] = 1.26482592e-18;
    veN2[28] = 1.31167132e-18;
    veN2[29] = 1.35851672e-18;
    veN2[30] = 1.40536213e-18;
    veN2[31] = 1.45220753e-18;
    veN2[32] = 1.49905294e-18;

    veO2[0] = 0.00000000e+00;
    veO2[1] = 3.13821409e-20; 
    veO2[2] = 6.27642817e-20;
    veO2[3] = 9.41464226e-20;
    veO2[4] = 1.25528563e-19;
    veO2[5] = 1.56910704e-19;
    veO2[6] = 1.88292845e-19;
    veO2[7] = 2.19674986e-19;
    veO2[8] = 2.51057127e-19;
    veO2[9] = 2.82439268e-19;
    veO2[10] = 3.13821409e-19; 
    veO2[11] = 3.45203549e-19;
    veO2[12] = 3.76585690e-19;
    veO2[13] = 4.07967831e-19;
    veO2[14] = 4.39349972e-19;
    veO2[15] = 4.70732113e-19;
    veO2[16] = 5.02114254e-19;
    veO2[17] = 5.33496395e-19;
    veO2[18] = 5.64878535e-19;
    veO2[19] = 5.96260676e-19;
    veO2[20] = 6.27642817e-19;
    veO2[21] = 6.59024958e-19;
    veO2[22] = 6.90407099e-19;
    veO2[23] = 7.21789240e-19;
    veO2[24] = 7.53171381e-19;
    veO2[25] = 7.84553521e-19;


    // Dissociation energy
    const double de[m_thermo.nMolecules()] = {1.0409017238088574e-18,
	1.5636156480913654e-18, 8.196091362099268e-19};
    for (int i=0; i<m_thermo.nMolecules(); ++i) {
	std::cout << "de: " << de[i] << std::endl;
    }

    // Partial number densities
    double nd[m_thermo.nSpecies()] = {};
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
        nd[i] = p_X[i] * n;
	std::cout << "nd: " << nd[i] << std::endl;
    }

    // Calculate the non-equilibrium partition functions
    double z_vibr_T_Tv_NO = 0.;
    for (int i=0; i<lNO; ++i) {

	// harmonic oscillator
        z_vibr_T_Tv_NO += exp(-veNO[i]/(KB*Tv[0]));

	// anharmonic oscillator
        //z_vibr_T_Tv_NO += exp(-(veNO[i]-i*veNO[0])/(KB*T)-
	//		                i*veNO[0]/(KB*Tv[0]));
    }
    double z_vibr_T_Tv_N2 = 0.;
    for (int i=0; i<lN2; ++i) {
        z_vibr_T_Tv_N2 += exp(-veN2[i]/(KB*Tv[1]));
        //z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/(KB*T)-
	//		                i*veN2[0]/(KB*Tv[1]));
    }
    double z_vibr_T_Tv_O2 = 0.;
    for (int i=0; i<lO2; ++i) {
        z_vibr_T_Tv_O2 += exp(-veO2[i]/(KB*Tv[2]));
        //z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/(KB*T)-
	//		                i*veO2[0]/(KB*Tv[2]));
    }
    std::cout << "zNO: " << z_vibr_T_Tv_NO << " " << 
	         "zN2: " << z_vibr_T_Tv_N2 << " " << 
		 "zO2: " << z_vibr_T_Tv_O2 << std::endl;


    // Compute the vibrational populations
    // according to Treanor-Marrone or Boltzmann
    double nNO[lNO];
    double nN2[lN2];
    double nO2[lO2];

    std::cout << " Vibrational population ... " << std::endl;
    for (int i=0; i<lNO; i++) {
        nNO[i] = nd[2] * exp(-veNO[i] / (KB * Tv[0])) / z_vibr_T_Tv_NO;
        //nNO[i] = nd[2] * exp(-(veNO[i] - i * veNO[0]) / KBT 
	//    - i * veNO[0] / (KB * Tv[0])) / z_vibr_T_Tv_NO;
        std::cout << i << " " << nNO[i] << std::endl;
    }
    for (int i=0; i<lN2; i++) {
        nN2[i] = nd[3] * exp(-veN2[i] / (KB * Tv[1])) / z_vibr_T_Tv_N2;
        //nN2[i] = nd[3] * exp(-(veN2[i] - i * veNO[0]) / KBT 
	//    - i * veN2[0] / (KB * Tv[1])) / z_vibr_T_Tv_N2;
        std::cout << i << " " << nN2[i] << std::endl;
    }
    for (int i=0; i<lO2; i++) {
        nO2[i] = nd[4] * exp(-veO2[i] / (KB * Tv[2])) / z_vibr_T_Tv_O2;
        //nO2[i] = nd[4] * exp(-(veO2[i] - i * veO2[0]) / KBT 
	//    - i * veO2[0] / (KB * Tv[2])) / z_vibr_T_Tv_O2;
        std::cout << i << " " << nO2[i] << std::endl;
    }

    // Masses 
    double m[5], d[5];
    m[0] = 2.3258672171567998e-26; // N
    m[1] = 2.6567628316576e-26;    // O
    m[2] = 4.9826300488143997e-26; // NO
    m[3] = 4.6517344343135997e-26; // N2
    m[4] = 5.3135256633152e-26;	   // O2

    // Diameters
    d[0] = 3.298e-10;
    d[1] = 2.75e-10;
    d[2] = 3.4061e-10;
    d[3] = 3.40385035355259e-10;
    d[4] = 3.5155e-10;
    
    // Compute state-to-state diss/rec coefficients, k_c,diss^d(0)
    // according to the Rigid-Sphere (RS) model.
    std::cout << " Dissociation STS rates ... " << std::endl;
    double kcd_no[5][lNO]; 
    double coll_mass = 0.;
    double diameter = 0.;
    double interNO [5][3] = {}; // 1st: interaction, 2nd: property
    interNO[0][0] = 1.82659281313541e-13; 	// A for N + NO 
    interNO[0][1] = 0.;			 	// n
    interNO[0][2] = 1.0423910520000002e-18;	// Ea
    interNO[1][0] = 1.82659281313541e-13; 	// A for O + NO 
    interNO[1][1] = 0.;			 	// n
    interNO[1][2] = 1.0423910520000002e-18;	// Ea
    interNO[2][0] = 1.82659281313541e-13; 	// A for NO + NO 
    interNO[2][1] = 0.;			 	// n
    interNO[2][2] = 1.0423910520000002e-18;	// Ea
    interNO[3][0] = 1.16237724472253e-08; 	// A for N2 + NO 
    interNO[3][1] = -1.6;			 	// n
    interNO[3][2] = 1.5628962528000002e-18;	// Ea
    interNO[4][0] = 2e-10;		 	// A for O2 + NO 
    interNO[4][1] = -1.;			 	// n
    interNO[4][2] = 1.043e-18;			// Ea
    double A, nn, Ea = 0.;
    for (int i=0; i<5; ++i) { // N O NO N2 O2
	coll_mass = (m[i] * m[2]) / (m[i] + m[2]);
	diameter = 0.5 * (d[i] + d[2]);
        for (int j=0; j<lNO; ++j) {

            // Treanor-Marrone, 3T model
	    A = interNO[i][0];
	    nn = interNO[i][1];
	    Ea = interNO[i][2];
	    kcd_no[i][j] = k_Arrhenius(T, A, nn, Ea) 
	        * Z_diss(T, 3.*T, veNO, j);

            // Treanor-Marrone, D/6K model
	    kcd_no[i][j] = k_Arrhenius(T, A, nn, Ea) 
	        * Z_diss(T, de[0]/(6.*KB), veNO, j);
	    
 	    // Rigid-sphere (RS) model	    
            kcd_no[i][j] = k_diss_RS(T, coll_mass, diameter, de[0], veNO[j],
	        /*center_of_mass=*/true);

 	    std::cout << i << " " << j << " " << T << " " << coll_mass 
		   << " " << diameter << " " << de[2] << " " << veNO[j] 
		   << " " << kcd_no[i][j] << std::endl;
	}
    }

    double interN2 [5][3] = {}; 
    interN2[0][0] = 2.657e-08; 			// A for N + N2
    interN2[0][1] = -1.6;			 
    interN2[0][2] = 1.5628962528000002e-18;
    interN2[1][0] = 4.98161676309657e-08; 	// A for O + N2 
    interN2[1][1] = -1.6;			 
    interN2[1][2] = 1.5628962528000002e-18;	
    interN2[2][0] = 1.16237724472253e-08; 	// A for NO + N2
    interN2[2][1] = -1.6;			 
    interN2[2][2] = 1.5628962528000002e-18;
    interN2[3][0] = 6.144e-09; 			// A for N2 + N2
    interN2[3][1] = -1.6;		
    interN2[3][2] = 1.5628962528000002e-18;
    interN2[4][0] = 1.16237724472253e-08; 	// A for O2 + N2
    interN2[4][1] = -1.6;		
    interN2[4][2] = 1.5628962528000002e-18;	
    double kcd_n2[5][lN2]; 
    for (int i=0; i<5; ++i) { // N O NO N2 O2
	coll_mass = (m[i] * m[3]) / (m[i] + m[3]);
	diameter = 0.5 * (d[i] + d[3]);
        for (int j=0; j<lN2; ++j) {

	    A = interN2[i][0];
	    nn = interN2[i][1];
	    Ea = interN2[i][2];
	    kcd_n2[i][j] = k_Arrhenius(T, A, nn, Ea) * Z_diss(T, 3*T, veN2, j);

	    kcd_n2[i][j] = k_Arrhenius(T, A, nn, Ea) 
	        * Z_diss(T, de[1]/(6.*KB), veN2, j);
	    
            kcd_n2[i][j] = k_diss_RS(T, coll_mass, diameter, de[1], veN2[j],
	        /*center_of_mass=*/true);

 	    std::cout << i << " " << j << " " << T << " " << coll_mass 
		   << " " << diameter << " " << de[1] << " " << veN2[j] 
		   << " " << kcd_n2[i][j] << std::endl;
	}
    }

    double interO2 [5][3] = {}; 
    interO2[0][0] = 3.32107784206438e-09;	// A for N + O2
    interO2[0][1] = -1.5;			 
    interO2[0][2] = 8.24938614e-19;
    interO2[1][0] = 1.66053892103219e-08; 	// A for O + O2 
    interO2[1][1] = -1.5;			 
    interO2[1][2] = 8.24938614e-19;	
    interO2[2][0] = 3.32107784206438e-09; 	// A for NO + O2
    interO2[2][1] = -1.5;			 
    interO2[2][2] = 8.24938614e-19;
    interO2[3][0] = 1.16237724472253e-08; 	// A for N2 + O2
    interO2[3][1] = -1.6;		
    interO2[3][2] = 1.5628962528000002e-18;
    interO2[4][0] = 3.32107784206438e-09; 	// A for O2 + O2
    interO2[4][1] = -1.5;		
    interO2[4][2] = 8.24938614e-19;	
    double kcd_o2[5][lO2]; 
    for (int i=0; i<5; ++i) { // N O NO N2 O2
	coll_mass = (m[i] * m[4]) / (m[i] + m[4]);
	diameter = 0.5 * (d[i] + d[4]);
        for (int j=0; j<lO2; ++j) {

	    A = interO2[i][0];
	    nn = interO2[i][1];
	    Ea = interO2[i][2];
	    kcd_o2[i][j] = k_Arrhenius(T, A, nn, Ea) 
                * Z_diss(T, 3.*T, veO2, j);
	
	    kcd_o2[i][j] = k_Arrhenius(T, A, nn, Ea) 
	        * Z_diss(T, de[2]/(6.*KB), veO2, j);

            kcd_o2[i][j] = k_diss_RS(T, coll_mass, diameter, de[2], veO2[j],
	        /*center_of_mass=*/true);

 	    std::cout << i << " " << j << " " << T << " " << coll_mass 
		   << " " << diameter << " " << de[2] << " " << veO2[j] 
		   << " " << kcd_o2[i][j] << std::endl;
	}
    }

    // Compute the rotational partition function
    // according to a simplified formula ... TODO: better formula
    double Be[3] = {199.8, 143.77, 167.20}; // m^-1
    int sigma[3] = {2, 2, 1};
    double theta_r[3] = {};
    double z_rot[3] = {};
    for (int i=0; i<3; ++i) {
        theta_r[i] = Be[i]*HP*C0/KB;
        z_rot[i] = T/(sigma[i]*theta_r[i]);
	std::cout << "theta & zrot:" 
	 	  << theta_r[i] << " " << z_rot[i] << std::endl;
    }

    // Compute STS recombination coefficients, k_rec,c^d(0)
    std::cout << " Recombination STS rates ... " << std::endl;
    double fac0 = HP*HP*HP * pow(2*PI*KB*T,-1.5);
    double fac1 = (m[2]/(m[0]*m[1])) * fac0 * z_rot[0];

    double kcr_no[5][lNO]; 
    for (int i=0; i<5; ++i) {
        for (int j=0; j<lNO; ++j) {
            kcr_no[i][j] = kcd_no[i][j] * exp(-(veNO[j]-de[0])/KBT) * fac1;
	    std::cout << i << " " << j << " " << kcr_no[i][j] << std::endl;
	}
    }

    fac1 = (m[3]/(m[0]*m[0])) * fac0 * z_rot[1];
    double kcr_n2[5][lN2]; 
    for (int i=0; i<5; ++i) {
        for (int j=0; j<lN2; ++j) {
            kcr_n2[i][j] = kcd_n2[i][j] * exp(-(veN2[j]-de[1])/KBT) * fac1;
	    std::cout << i << " " << j << " " << kcr_n2[i][j] << std::endl;
	}
    }

    fac1 = (m[4]/(m[1]*m[1])) * fac0 * z_rot[2];
    double kcr_o2[5][lO2]; 
    for (int i=0; i<5; ++i) {
        for (int j=0; j<lO2; ++j) {
            kcr_o2[i][j] = kcd_o2[i][j] * exp(-(veO2[j]-de[2])/KBT) * fac1;
	    std::cout << i << " " << j << " " << kcr_o2[i][j] << std::endl;
	}
    }

    // Compute MT diss/rec coefficients, k_c^d(0), eq. 3.94.
    double kd_no[m_thermo.nSpecies()];
    double kr_no[m_thermo.nSpecies()];
    double kd_n2[m_thermo.nSpecies()];
    double kr_n2[m_thermo.nSpecies()];
    double kd_o2[m_thermo.nSpecies()];
    double kr_o2[m_thermo.nSpecies()];

    // Compute the summation term in eq. 3.94.
    // NO
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
	fac0 = 0., fac1 = 0.; 
        for (int j=0; j<lNO; ++j) {

	    // Using the Boltzmann distribution
	    fac0 += nNO[j] * kcd_no[i][j];

	    // Using the generalized Treanor-Marrone distribution	
            //fac0 += exp(-((veNO[j]-j*veNO[0])/KBT) -
            //                       j*veNO[0]/(KB*Tv[0])) * kcd_no[i][j];
	    fac1 += kcr_no[i][j];
        }
	kd_no[i] = fac0 / nd[2];
        //kd_no[i] = fac0 / z_vibr_T_Tv_NO;
        kr_no[i] = fac1;
	std::cout << i << " kdNO: " << kd_no[i] 
		       << " krNO: " << kr_no[i] << std::endl;
    }

    // N2
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
	fac0 = 0., fac1 = 0.; 
        for (int j=0; j<lN2; ++j) {

	    fac0 += nN2[j] * kcd_n2[i][j];

            //fac0 += exp(-((veN2[j]-j*veN2[0])/KBT) -
            //                       j*veN2[0]/(KB*Tv[1])) * kcd_n2[i][j];
	    fac1 += kcr_n2[i][j];
        }
	kd_n2[i] = fac0 / nd[3];
        //kd_n2[i] = fac0 / z_vibr_T_Tv_N2;
        kr_n2[i] = fac1;
	std::cout << i << " kdN2: " << kd_n2[i] 
		       << " krN2: " << kr_n2[i] << std::endl;
    }

    // O2
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
	fac0 = 0., fac1 = 0.; 
        for (int j=0; j<lO2; ++j) {

	    fac0 += nO2[j] * kcd_o2[i][j];

            //fac0 += exp(-((veO2[j]-j*veO2[0])/KBT) -
            //                       j*veO2[0]/(KB*Tv[2])) * kcd_o2[i][j];
	    fac1 += kcr_o2[i][j];
        }
	kd_o2[i] = fac0 / nd[4];
        //kd_o2[i] = fac0 / z_vibr_T_Tv_O2;
        kr_o2[i] = fac1;
	std::cout << i << " kdO2: " << kd_o2[i] 
		       << " krO2: " << kr_o2[i] << std::endl;
    }

    // N O NO N2 O2
    double R_react_23[m_thermo.nMolecules()] = {};

    // NO
    R_react_23[0] = nd[0] * (nd[0]*nd[1]*kr_no[0]-nd[2]*kd_no[0]) +
                    nd[1] * (nd[0]*nd[1]*kr_no[1]-nd[2]*kd_no[1]) +
                    nd[2] * (nd[0]*nd[1]*kr_no[2]-nd[2]*kd_no[2]) +
                    nd[3] * (nd[0]*nd[1]*kr_no[3]-nd[2]*kd_no[3]) +
                    nd[4] * (nd[0]*nd[1]*kr_no[4]-nd[2]*kd_no[4]) ;
    // N2
    R_react_23[1] = nd[0] * (nd[0]*nd[0]*kr_n2[0]-nd[3]*kd_n2[0]) +
                    nd[1] * (nd[0]*nd[0]*kr_n2[1]-nd[3]*kd_n2[1]) +
                    nd[2] * (nd[0]*nd[0]*kr_n2[2]-nd[3]*kd_n2[2]) +
                    nd[3] * (nd[0]*nd[0]*kr_n2[3]-nd[3]*kd_n2[3]) +
                    nd[4] * (nd[0]*nd[0]*kr_n2[4]-nd[3]*kd_n2[4]) ;
    // O2
    R_react_23[2] = nd[0] * (nd[1]*nd[1]*kr_o2[0]-nd[4]*kd_o2[0]) +
                    nd[1] * (nd[1]*nd[1]*kr_o2[1]-nd[4]*kd_o2[1]) +
                    nd[2] * (nd[1]*nd[1]*kr_o2[2]-nd[4]*kd_o2[2]) +
                    nd[3] * (nd[1]*nd[1]*kr_o2[3]-nd[4]*kd_o2[3]) +
                    nd[4] * (nd[1]*nd[1]*kr_o2[4]-nd[4]*kd_o2[4]) ;

    // Compute species concentrations (mol/m^3)
    Map<ArrayXd>(p_r_23, m_thermo.nMolecules()) = 0.;

    // TODO: check dimensions!
    for (int i=0; i<m_thermo.nMolecules(); ++i) {
        p_r_23[i] = R_react_23[i];
	std::cout <<  " R_23["<<i<<"]" << p_r_23[i] << std::endl;
    }
}

//==============================================================================

    double Kinetics::k_diss_RS(double T, double coll_mass, double diameter, 
        double diss_energy, double vibr_energy, bool center_of_mass) 
    {
        return 8. * integral_diss_RS(T, 0, coll_mass, diameter, diss_energy, 
	    vibr_energy, center_of_mass);
    }

//==============================================================================

    double Kinetics::k_diss_VSS(double T, double coll_mass, double vss_c_cs, 
        double vss_omega, double diss_energy, double vibr_energy, 
	    bool center_of_mass) 
    {
        return 8. * integral_diss_VSS(T, 0, coll_mass, vss_c_cs, vss_omega, 
	    diss_energy, vibr_energy, center_of_mass);
  }

//==============================================================================

    double Kinetics::integral_diss_RS(double T, int degree,
        double coll_mass, double diameter, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {
        double conversion = sqrt(2. * KB * T / coll_mass);

        auto integrand = [=](double g)
        {
	    return pow(g, 2 * degree + 3) * 
	        crosssection_diss_RS(conversion * g, coll_mass, diameter, 
		    diss_energy, vibr_energy, center_of_mass) * exp(-g * g); 
	};

        return sqrt(KB * T / (2 * PI * coll_mass)) * 
	    Numerics::integrate_semi_inf(integrand, min_vel_diss(coll_mass, 
                diss_energy, vibr_energy) / conversion);
    }

//==============================================================================

    double Kinetics::integral_diss_VSS(double T, int degree, double coll_mass, 
        double vss_c_cs, double vss_omega, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {

        double conversion = sqrt(2 * KB * T / coll_mass);

        auto integrand = [=](double g)
	{
            return pow(g, 2 * degree + 3) * crosssection_diss_VSS(
	        conversion * g, coll_mass, vss_c_cs, vss_omega, 
	            diss_energy, vibr_energy, center_of_mass) * 
		        exp(-g * g); 
	};

        return sqrt(KB * T / (2. * PI * coll_mass)) * 
	    integrate_semi_inf(integrand, min_vel_diss(coll_mass, 
	        diss_energy, vibr_energy) / conversion);
    }

//==============================================================================

    double Kinetics::crosssection_diss_RS(double rel_vel, double coll_mass, 
        double diameter, double diss_energy, double vibr_energy, 
	    bool center_of_mass) 
    {
        return crosssection_elastic_RS(diameter) * 
	    p_probability_diss(rel_vel, coll_mass, diss_energy, vibr_energy, 
	        center_of_mass);
    }

//==============================================================================

    double Kinetics::crosssection_diss_VSS(double rel_vel, double coll_mass, 
        double vss_c_cs, double vss_omega, double diss_energy, 
	    double vibr_energy, bool center_of_mass) 
    {   
        return crosssection_elastic_VSS(rel_vel, vss_c_cs, vss_omega) *
	    p_probability_diss(rel_vel, coll_mass, diss_energy, vibr_energy, 
	        center_of_mass);
    }

//==============================================================================

    double Kinetics::p_probability_diss(double rel_vel, double coll_mass, 
        double diss_energy, double vibr_energy, bool center_of_mass) 
    {

        double energy = vibr_energy + rel_vel * rel_vel * coll_mass / 2.;

        if (energy < diss_energy) {
            return 0.0;
        } else if (center_of_mass) {
            return 1.0 - diss_energy / energy;
        } else {
            return 1.0;
        }
    }

//==============================================================================

    double Kinetics::min_vel_diss(double coll_mass, double diss_energy, 
        double vibr_energy) 
    {
        return sqrt(2. * (diss_energy - vibr_energy) / coll_mass);
    }

//==============================================================================

    double Kinetics::crosssection_elastic_RS(double diameter) 
    {
        return PI * diameter * diameter;
    }

//==============================================================================

    double Kinetics::crosssection_elastic_VSS(double rel_vel, double vss_c_cs, 
	double vss_omega) 
    {
        return vss_c_cs * pow(rel_vel, 1. - 2. * vss_omega);
    }

//==============================================================================

    double Kinetics::crosssection_elastic_VSS(double rel_vel, 
	double coll_mass, double vss_c, double vss_omega) 
    {
        return vss_c * pow(coll_mass * rel_vel * rel_vel / (2. * KB), 
	    -vss_omega);
    }

//==============================================================================

    double Kinetics::k_Arrhenius(double T, double A, double n, double Ea) 
    {

        std::cout << " In k_Arrhenius ... " << std::endl;    
        std::cout << T << std::endl;    
        std::cout << A << std::endl;    
        std::cout << n << std::endl;    
        std::cout << Ea << std::endl;    
        return A * pow(T, n) * exp(-Ea / (KB * T));
    }

//==============================================================================

   // Compute the non-equilibrium factor using the Treanor-Marrone model
   double Kinetics::Z_diss(double T, double U, const double *ve, int i) 
   {
       return p_Z_vibr_eq(T, ve) * exp(ve[i] * (1./T + 1./U) / KB) 
	       / p_Z_vibr_eq(-U, ve);
   }

//==============================================================================

  // Computes the equilibrium vibrational partition function
  double Kinetics::p_Z_vibr_eq(double T, const double *ve) 
  {
    double p_Z_vibr_eq = 0.;
    int size = sizeof(ve)/sizeof(ve[0]);
    for (int i=0; i<size; ++i) 
        p_Z_vibr_eq += exp(-ve[i] / (KB * T));
        
    return p_Z_vibr_eq;
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
    
    for(int i=0; i < mp_rates->irrReactions().size(); ++i)
        mp_ropb[mp_rates->irrReactions()[i]] = 0.0;
    
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
