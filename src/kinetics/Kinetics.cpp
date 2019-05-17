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
      //mp_r_vt(NULL),
      //mp_r_vv(NULL),
      //mp_r_diss_rec(NULL)
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
    //if (mp_r_vt != NULL)
    //    delete [] mp_r_vt;
    //if (mp_r_vv != NULL)
    //    delete [] mp_r_vv;
    //if (mp_r_diss_rec != NULL)
    //    delete [] mp_r_diss_rec;
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

    mp_r          = new double [m_thermo.nSpecies()];
    //mp_r_vt       = new double [m_thermo.nSpecies()];
    //mp_r_vv       = new double [m_thermo.nSpecies()];
    //mp_r_diss_rec = new double [m_thermo.nSpecies()];
    mp_r_23       = new double [m_thermo.nSpecies()];
    
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
// The whole block of functions has been moved in R.h and R.cpp 

void Kinetics::R(double* const p_r)
{
    // Sum all contributions from every process
    std::fill(p_r, p_r+m_thermo.nSpecies(), 0.0);
    R_22(mp_r_23);
    R_23(mp_r_23);

    // Multiply by species molecular weights
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_r[i] *= m_thermo.speciesMw(i);
}

//==============================================================================

// Compute the production term R_22 due to exchange reactions
// The dimensions of the source production term, p_r_diss_rec, are
// given by the number of species present in the mixture.
// For the moment, I just do this for air5, MT.
void Kinetics::R_22(double* const p_r_22)
{
    // Retrieve all necessary temperatures
    double T = m_thermo.T();
    double Tv[3];
    for (int i=0; i<3; ++i)
        Tv[i] = m_thermo.Tvs(i); // NO, N2, O2
    Tv[0] = 4000;
    Tv[1] = 4000;
    Tv[2] = 4000;
    std::cout << T     << " "
	      << Tv[0] << " "  /*NO*/
	      << Tv[1] << " "  /*N2*/
	      << Tv[2] << std::endl; /*O2*/

    // The following data section is brutally hard-coded 
    // and should be put somewhere else
    // as well as being flexible and general TODO

    // Number of vibrational levels
    int lN2 = 48;
    int lO2 = 36;
    int lNO = 38;

    // Vibrational energy arrays
    static double ve[122];
    static double veN2[48];
    static double veO2[36];
    static double veNO[38];

    ve[1] = 6.96372e-20;
    ve[2] = 1.1535e-19 ;
    ve[3] = 1.60493e-19;
    ve[4] = 2.05065e-19;
    ve[5] = 2.49065e-19;
    ve[6] = 2.92494e-19;
    ve[7] = 3.35349e-19;
    ve[8] = 3.77629e-19;
    ve[9] = 4.19334e-19;
    ve[10] = 4.60463e-19;
    ve[11] = 5.01013e-19;
    ve[12] = 5.40983e-19;
    ve[13] = 5.80372e-19;
    ve[14] = 6.19178e-19;
    ve[15] = 6.57399e-19;
    ve[16] = 6.95033e-19;
    ve[17] = 7.32077e-19;
    ve[18] = 7.68531e-19;
    ve[19] = 8.0439e-19;
    ve[20] = 8.39654e-19;
    ve[21] = 8.74319e-19;
    ve[22] = 9.08383e-19;
    ve[23] = 9.41842e-19;
    ve[24] = 9.74695e-19;
    ve[25] = 1.00694e-18;
    ve[26] = 1.03857e-18;
    ve[27] = 1.06958e-18;
    ve[28] = 1.09997e-18;
    ve[29] = 1.12974e-18;
    ve[30] = 1.15889e-18;
    ve[31] = 1.1874e-18;
    ve[32] = 1.21528e-18;
    ve[33] = 1.24252e-18;
    ve[34] = 1.26911e-18;
    ve[35] = 1.29507e-18;
    ve[36] = 1.32037e-18;
    ve[37] = 1.34501e-18;
    ve[38] = 1.369e-18  ;
    ve[39] = 1.39232e-18;
    ve[40] = 1.41497e-18;
    ve[41] = 1.43695e-18;
    ve[42] = 1.45825e-18;
    ve[43] = 1.47887e-18;
    ve[44] = 1.49879e-18;
    ve[45] = 1.51803e-18;
    ve[46] = 1.53656e-18;
    ve[47] = 1.55439e-18;
    ve[48] = 1.56354e-20;
    ve[49] = 4.6552e-20;
    ve[50] = 7.70004e-20;
    ve[51] = 1.06985e-19;
    ve[52] = 1.3651e-19;
    ve[53] = 1.65578e-19;
    ve[54] = 1.94192e-19;
    ve[55] = 2.22354e-19;
    ve[56] = 2.50065e-19;
    ve[57] = 2.77327e-19;
    ve[58] = 3.04138e-19;
    ve[59] = 3.305e-19  ;
    ve[60] = 3.56411e-19;
    ve[61] = 3.81869e-19;
    ve[62] = 4.06872e-19;
    ve[63] = 4.31417e-19;
    ve[64] = 4.55501e-19;
    ve[65] = 4.7912e-19;
    ve[66] = 5.02269e-19;
    ve[67] = 5.24943e-19;
    ve[68] = 5.47135e-19;
    ve[69] = 5.6884e-19;
    ve[70] = 5.90051e-19;
    ve[71] = 6.10759e-19;
    ve[72] = 6.30957e-19;
    ve[73] = 6.50635e-19;
    ve[74] = 6.69784e-19;
    ve[75] = 6.88393e-19;
    ve[76] = 7.06453e-19;
    ve[77] = 7.23952e-19;
    ve[78] = 7.40877e-19;
    ve[79] = 7.57217e-19;
    ve[80] = 7.72958e-19;
    ve[81] = 7.88086e-19;
    ve[82] = 8.02588e-19;
    ve[83] = 8.16447e-19;
    ve[84] = 1.88431e-20;
    ve[85] = 5.61105e-20;
    ve[86] = 9.28207e-20;
    ve[87] = 1.28975e-19;
    ve[88] = 1.64575e-19;
    ve[89] = 1.99621e-19;
    ve[90] = 2.34116e-19;
    ve[91] = 2.68059e-19;
    ve[92] = 3.01454e-19;
    ve[93] = 3.343e-19  ;
    ve[94] = 3.666e-19  ;
    ve[95] = 3.98354e-19;
    ve[96] = 4.29564e-19;
    ve[97] = 4.60232e-19;
    ve[98] = 4.90357e-19;
    ve[99] = 5.19943e-19;
    ve[100] = 5.4899e-19;
    ve[101] = 5.77499e-19;
    ve[102] = 6.05472e-19;
    ve[103] = 6.3291e-19 ;
    ve[104] = 6.59815e-19;
    ve[105] = 6.86187e-19;
    ve[106] = 7.12028e-19;
    ve[107] = 7.3734e-19 ;
    ve[108] = 7.62123e-19;
    ve[109] = 7.86379e-19;
    ve[110] = 8.10109e-19;
    ve[111] = 8.33315e-19;
    ve[112] = 8.55998e-19;
    ve[113] = 8.78159e-19;
    ve[114] = 8.99799e-19;
    ve[115] = 9.2092e-19 ;
    ve[116] = 9.41523e-19;
    ve[117] = 9.6161e-19 ;
    ve[118] = 9.81182e-19;
    ve[119] = 1.00024e-18;
    ve[120] = 1.01878e-18;
    ve[121] = 1.03682e-18;

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

    // Dissociation energy (ground level)
    const double deNO = 1.0409017238088574e-18; //[J]
    const double deN2 = 1.5636156480913654e-18; //[J]
    const double deO2 = 8.196091362099268e-19;  //[J]

    // Special case of no reactions
    if (nReactions() == 0) {
        std::fill(p_r_22, p_r_22 + m_thermo.nSpecies(), 0);
        return;
    }

    // If there are no diss/rec reactions ... TODO

    // Compute species concentrations (mol/m^3)
    Map<ArrayXd>(p_r_22, m_thermo.nSpecies()) = 0.;
    //    (m_thermo.numberDensity() / NA) *
    //    Map<const ArrayXd>(m_thermo.X(), m_thermo.nSpecies());

    const double * p_Y = m_thermo.Y();
    const double * p_X = m_thermo.X();
    const double   n   = m_thermo.numberDensity();
    const double   rho = m_thermo.density();
    std::cout << "Mass and Molar fractions" << std::endl;
    std::cout << p_Y[0] << " " << p_Y[1] << " " << p_Y[2] << " " << p_Y[3] 
	      << p_Y[4] << "\n"  
	      << p_X[0] << " " << p_X[1] << " " << p_X[2] << " " << p_X[3] 
              << p_X[4] << "\n"	
	      << n   << " " 
	      << rho << std::endl;

    // Compute partial number densities,
    // the order is given by input file, air_5.xml
    double nd_n  = p_X[0] / n;
    double nd_o  = p_X[1] / n;
    double nd_no = p_X[2] / n;
    double nd_n2 = p_X[3] / n;
    double nd_o2 = p_X[4] / n;
    std::cout << "Number densities" << std::endl;
    std::cout << nd_n  << " " 
	      << nd_o  << " " 	
	      << nd_no << " " 
	      << nd_n2 << " " 
	      << nd_o2 << std::endl;

    // Calculate the non-equilibrium partition functions
    double z_vibr_T_Tv_NO = 0.;
    for (int i=0; i<lNO; ++i) {
        z_vibr_T_Tv_NO += exp(-(veNO[i]-i*veNO[0])/(KB*T)-
			                i*veNO[0]/(KB*Tv[0]));
    }
    double z_vibr_T_Tv_N2 = 0.;
    for (int i=0; i<lN2; ++i) {
        z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/(KB*T)-
			                i*veN2[0]/(KB*Tv[1]));
    }
    double z_vibr_T_Tv_O2 = 0.;
    for (int i=0; i<lO2; ++i) {
        z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/(KB*T)-
			                i*veO2[0]/(KB*Tv[2]));
    }
    std::cout << "Non-eq. partition functions" << std::endl;
    std::cout << z_vibr_T_Tv_NO << " " 
	      << z_vibr_T_Tv_N2 << " " 	
	      << z_vibr_T_Tv_O2 << std::endl;

    // Compute state-to-state diss/rec coefficients, k_c,diss^d(0)
    // according to the simplest Rigid Sphere (RS) model
    // TODO: add better models!

    // Mass and diameter
    double m_NO = 4.9826300488143997e-26; 	// kg
    double d_NO = 3.4061e-10; 			// m
    double m_N2 = 4.6517344343135997e-26;
    double d_N2 = 3.40385035355259e-10; 
    double m_O2 = 5.3135256633152e-26;
    double d_O2 = 3.5155e-10; 
    double m_N  = 2.3258672171567998e-26;
    double m_O  = 2.6567628316576e-26;

    double m[5], d[5];
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

    // TODO: Compute STS forward rates ...
    // not clear how to do it!
    double k_N2_NO[lO2][lNO];
    double k_O2_NO[lO2][lNO];
    double k_NO_N2[lO2][lNO];
    double k_NO_O2[lO2][lNO];

    // TODO: here is the last and most delicate part:
    // calculation of STS forward rates according to some model!
    double k_N2_O[lN2];
    // TODO: this should be selected in the input file, so it 
    // should be made polymorphic ...
    // TODO: here the Arrhenius formula has been hard-coded
    // but should be called from somewhere else!
    std::string models_k_exch = "arrh_park";
    for (int i=0; i<lN2; ++i) { // so it does not depend on i?
	//k_N2_O[i] = k_exch(T, i, models_k_exch);
	k_N2_O[i] = 1.0627449094606e-12/*8e-17*/ 
		    * pow(T, -1/*0.*/) 
		    * exp(-5.297555584800001e-19/*5.175e-19*/ / (KB * T));
    }
    double k_O2_N[lO2];
    for (int i=0; i<lO2; ++i) {
	//k_O2_N[i] = k_exch(T, i, models_k_exch);
	k_O2_N[i] = 4e-15 
		    * pow(T, -0.39) 
		    * exp(-2e-20 / (KB * T));
    }

    // Compute the rotational partition function
    // according to a simplified formula ... TODO: better formula
    double Be[3] = {1.998, 1.4377, 1.6720}; // *100 // m^-1
    int sigma[3] = {2, 2, 1};
    double theta_r[3] = {0, 0, 0};
    double z_rot[3] = {0, 0, 0};
    for (int i=0; i<3; ++i) {
        theta_r[i] = Be[i]*100*HP*C0/KB;
        z_rot[i] = T/(sigma[i]*theta_r[i]);
    }

    // Compute the equilibrium factors
    // TODO: check order of elementy
    double exp_no_n2[lN2][lNO] = {}; 
    for (int i=0; i<lN2; ++i) {
        for (int j=0; j<lNO; ++j) {
	    exp_no_n2[i][j] = exp(veNO[j] - veN2[i]);
	}
    }
    double exp_no_o2[lO2][lNO] = {}; 
    for (int i=0; i<lO2; ++i) {
        for (int j=0; j<lNO; ++j) {
	    exp_no_o2[i][j] = exp(veNO[j] - veO2[i]);
	}
    }

    double Kz_n2[lN2][lNO] = {}; 
    for (int i=0; i<lN2; ++i) {
        for (int j=0; j<lNO; ++j) {
            Kz_n2[i][j] = pow((m_N2*m_O/(m_NO*m_N)),1.5) * z_rot[1]/z_rot[0] *
                exp_no_n2[i][j] /(KB*T) * exp((deN2-deNO)/KB*T);
	}
    }
    double Kz_o2[lO2][lNO] = {}; 
    for (int i=0; i<lO2; ++i) {
        for (int j=0; j<lNO; ++j) {
            Kz_o2[i][j] = pow((m_O2*m_N/(m_NO*m_O)),1.5) * z_rot[2]/z_rot[0] *
                exp_no_o2[i][j] /(KB*T) * exp((deO2-deNO)/KB*T);
	}
    }

    // Compute MT forward exchange coefficients, eq. 3.93.
    // N2
    double fac0f, fac0b = 0.;
    for (int i=0; i<lN2; ++i) {
        for (int j=0; j<lNO; ++j) {
	    // MT forward rates
            fac0f += exp(-((veN2[i]-i*veN2[0])/(KB*T)) -
	                            i*veN2[0]/(KB*Tv[1])) * k_N2_NO[i][j];
	    // STS backward rates
	    k_NO_N2[i][j] = Kz_n2[i][j] * k_N2_NO[i][j];
            fac0b += exp(-((veN2[i]-i*veN2[0])/(KB*T)) -
	                            i*veN2[0]/(KB*Tv[1])) * k_NO_N2[i][j];
	}
    }
    double k_exch_N2_NO = fac0f / z_vibr_T_Tv_N2;
    double k_exch_NO_N2 = fac0b / z_vibr_T_Tv_N2;

    // O2
    for (int i=0; i<lO2; ++i) {
	fac0f = fac0b = 0.;
        for (int j=0; j<lNO; ++j) {
	    // MT forward rates
            fac0f += exp(-((veO2[i]-i*veO2[0])/(KB*T)) -
                                   i*veO2[0]/(KB*Tv[2])) * k_O2_NO[i][j];
	    // STS backward rates
	    k_NO_O2[i][j] = Kz_o2[i][j] * k_NO_O2[i][j];
            fac0b += exp(-((veO2[i]-i*veO2[0])/(KB*T)) -
                                    i*veO2[0]/(KB*Tv[2])) * k_NO_O2[i][j];
        }
    }
    double k_exch_O2_NO = fac0f / z_vibr_T_Tv_O2;
    double k_exch_NO_O2 = fac0b / z_vibr_T_Tv_O2;

    // Here, as elsewhere, a certain fixed order of species is assumed.
    // Such order is given by the list of species in the input file:
    // N O NO N2 O2
    double R_exch_22[m_thermo.nSpecies()] = {0.,0.,0.,0.,0.};
    double nd[m_thermo.nSpecies()] = {nd_n, nd_o, nd_no, nd_n2, nd_o2};

    // N
    R_exch_22[0] = 0.;
    // O
    R_exch_22[1] = 0.;
    // NO
    R_exch_22[2] = (nd_n2 * nd_o * k_exch_N2_NO - nd_no * nd_n * k_exch_NO_N2) + 
	           (nd_o2 * nd_n * k_exch_O2_NO - nd_no * nd_o * k_exch_NO_O2);
    // N2
    R_exch_22[3] = nd_no * nd_n * k_exch_NO_N2 - nd_n2 * nd_o * k_exch_N2_NO;
    // O2
    R_exch_22[4] = nd_no * nd_o * k_exch_NO_O2 - nd_o2 * nd_n * k_exch_O2_NO;

    // TODO: check dimensions!
    for (int i = 0; i < m_thermo.nSpecies(); ++i) {
        p_r_22[i] = R_exch_22[i];
	std::cout << p_r_22[i] << std::endl;
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
      return 4 * PI * n * diameter*diameter * sqrt(KB * T / (2*PI * coll_mass));
  }

//==============================================================================

// Compute the production term R_23 due to dissociation/recombination
// The dimensions of the source production term, p_r_diss_rec, are
// given by the number of species present in the mixture.
// For the moment, I just do this for air5, MT.
void Kinetics::R_23(double* const p_r_23)
{
    // Retrieve all necessary temperatures
    double T = m_thermo.T();
    double Tv[3];
    for (int i=0; i<3; ++i)
        Tv[i] = m_thermo.Tvs(i); // NO, N2, O2
    Tv[0] = 4000;
    Tv[1] = 4000;
    Tv[2] = 4000;
    std::cout << T     << " "
	      << Tv[0] << " "  /*NO*/
	      << Tv[1] << " "  /*N2*/
	      << Tv[2] << std::endl; /*O2*/

    // The following data section is brutally hardcoded 
    // and clearly should be retrieved from Particle
    // as well as being flexible and general TODO

    // Number of vibrational levels
    int lN2 = 48;
    int lO2 = 36;
    int lNO = 38;

    // Vibrational energy arrays
    static double ve[122];
    static double veN2[48];
    static double veO2[36];
    static double veNO[38];

    ve[1] = 6.96372e-20;
    ve[2] = 1.1535e-19 ;
    ve[3] = 1.60493e-19;
    ve[4] = 2.05065e-19;
    ve[5] = 2.49065e-19;
    ve[6] = 2.92494e-19;
    ve[7] = 3.35349e-19;
    ve[8] = 3.77629e-19;
    ve[9] = 4.19334e-19;
    ve[10] = 4.60463e-19;
    ve[11] = 5.01013e-19;
    ve[12] = 5.40983e-19;
    ve[13] = 5.80372e-19;
    ve[14] = 6.19178e-19;
    ve[15] = 6.57399e-19;
    ve[16] = 6.95033e-19;
    ve[17] = 7.32077e-19;
    ve[18] = 7.68531e-19;
    ve[19] = 8.0439e-19;
    ve[20] = 8.39654e-19;
    ve[21] = 8.74319e-19;
    ve[22] = 9.08383e-19;
    ve[23] = 9.41842e-19;
    ve[24] = 9.74695e-19;
    ve[25] = 1.00694e-18;
    ve[26] = 1.03857e-18;
    ve[27] = 1.06958e-18;
    ve[28] = 1.09997e-18;
    ve[29] = 1.12974e-18;
    ve[30] = 1.15889e-18;
    ve[31] = 1.1874e-18;
    ve[32] = 1.21528e-18;
    ve[33] = 1.24252e-18;
    ve[34] = 1.26911e-18;
    ve[35] = 1.29507e-18;
    ve[36] = 1.32037e-18;
    ve[37] = 1.34501e-18;
    ve[38] = 1.369e-18  ;
    ve[39] = 1.39232e-18;
    ve[40] = 1.41497e-18;
    ve[41] = 1.43695e-18;
    ve[42] = 1.45825e-18;
    ve[43] = 1.47887e-18;
    ve[44] = 1.49879e-18;
    ve[45] = 1.51803e-18;
    ve[46] = 1.53656e-18;
    ve[47] = 1.55439e-18;
    ve[48] = 1.56354e-20;
    ve[49] = 4.6552e-20;
    ve[50] = 7.70004e-20;
    ve[51] = 1.06985e-19;
    ve[52] = 1.3651e-19;
    ve[53] = 1.65578e-19;
    ve[54] = 1.94192e-19;
    ve[55] = 2.22354e-19;
    ve[56] = 2.50065e-19;
    ve[57] = 2.77327e-19;
    ve[58] = 3.04138e-19;
    ve[59] = 3.305e-19  ;
    ve[60] = 3.56411e-19;
    ve[61] = 3.81869e-19;
    ve[62] = 4.06872e-19;
    ve[63] = 4.31417e-19;
    ve[64] = 4.55501e-19;
    ve[65] = 4.7912e-19;
    ve[66] = 5.02269e-19;
    ve[67] = 5.24943e-19;
    ve[68] = 5.47135e-19;
    ve[69] = 5.6884e-19;
    ve[70] = 5.90051e-19;
    ve[71] = 6.10759e-19;
    ve[72] = 6.30957e-19;
    ve[73] = 6.50635e-19;
    ve[74] = 6.69784e-19;
    ve[75] = 6.88393e-19;
    ve[76] = 7.06453e-19;
    ve[77] = 7.23952e-19;
    ve[78] = 7.40877e-19;
    ve[79] = 7.57217e-19;
    ve[80] = 7.72958e-19;
    ve[81] = 7.88086e-19;
    ve[82] = 8.02588e-19;
    ve[83] = 8.16447e-19;
    ve[84] = 1.88431e-20;
    ve[85] = 5.61105e-20;
    ve[86] = 9.28207e-20;
    ve[87] = 1.28975e-19;
    ve[88] = 1.64575e-19;
    ve[89] = 1.99621e-19;
    ve[90] = 2.34116e-19;
    ve[91] = 2.68059e-19;
    ve[92] = 3.01454e-19;
    ve[93] = 3.343e-19  ;
    ve[94] = 3.666e-19  ;
    ve[95] = 3.98354e-19;
    ve[96] = 4.29564e-19;
    ve[97] = 4.60232e-19;
    ve[98] = 4.90357e-19;
    ve[99] = 5.19943e-19;
    ve[100] = 5.4899e-19;
    ve[101] = 5.77499e-19;
    ve[102] = 6.05472e-19;
    ve[103] = 6.3291e-19 ;
    ve[104] = 6.59815e-19;
    ve[105] = 6.86187e-19;
    ve[106] = 7.12028e-19;
    ve[107] = 7.3734e-19 ;
    ve[108] = 7.62123e-19;
    ve[109] = 7.86379e-19;
    ve[110] = 8.10109e-19;
    ve[111] = 8.33315e-19;
    ve[112] = 8.55998e-19;
    ve[113] = 8.78159e-19;
    ve[114] = 8.99799e-19;
    ve[115] = 9.2092e-19 ;
    ve[116] = 9.41523e-19;
    ve[117] = 9.6161e-19 ;
    ve[118] = 9.81182e-19;
    ve[119] = 1.00024e-18;
    ve[120] = 1.01878e-18;
    ve[121] = 1.03682e-18;

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

    // Dissociation energy
    const double deNO = 1.0409017238088574e-18; //[J]
    const double deN2 = 1.5636156480913654e-18; //[J]
    const double deO2 = 8.196091362099268e-19;  //[J]

    // Special case of no reactions
    if (nReactions() == 0) {
        std::fill(p_r_23, p_r_23 + m_thermo.nSpecies(), 0);
        return;
    }

    // If there are no diss/rec reactions ... TODO

    // Compute species concentrations (mol/m^3)
    Map<ArrayXd>(p_r_23, m_thermo.nSpecies()) = 0.;
    //    (m_thermo.numberDensity() / NA) *
    //    Map<const ArrayXd>(m_thermo.X(), m_thermo.nSpecies());

    const double * p_Y = m_thermo.Y();
    const double * p_X = m_thermo.X();
    const double   n   = m_thermo.numberDensity();
    const double   rho = m_thermo.density();
    std::cout << "Mass and Molar fractions" << std::endl;
    std::cout << p_Y[0] << " " << p_Y[1] << " " << p_Y[2] << " " << p_Y[3] 
	      << p_Y[4] << "\n"  
	      << p_X[0] << " " << p_X[1] << " " << p_X[2] << " " << p_X[3] 
              << p_X[4] << "\n"	
	      << n   << " " 
	      << rho << std::endl;

    // Compute partial number densities,
    // the order is given by input file, air_5.xml
    double nd_n  = p_X[0] / n;
    double nd_o  = p_X[1] / n;
    double nd_no = p_X[2] / n;
    double nd_n2 = p_X[3] / n;
    double nd_o2 = p_X[4] / n;
    std::cout << "Number densities" << std::endl;
    std::cout << nd_n  << " " 
	      << nd_o  << " " 	
	      << nd_no << " " 
	      << nd_n2 << " " 
	      << nd_o2 << std::endl;

    // Calculate the non-equilibrium partition functions
    double z_vibr_T_Tv_NO = 0.;
    for (int i=0; i<lNO; ++i) {
        z_vibr_T_Tv_NO += exp(-(veNO[i]-i*veNO[0])/(KB*T)-
			                i*veNO[0]/(KB*Tv[0]));
    }
    double z_vibr_T_Tv_N2 = 0.;
    for (int i=0; i<lN2; ++i) {
        z_vibr_T_Tv_N2 += exp(-(veN2[i]-i*veN2[0])/(KB*T)-
			                i*veN2[0]/(KB*Tv[1]));
    }
    double z_vibr_T_Tv_O2 = 0.;
    for (int i=0; i<lO2; ++i) {
        z_vibr_T_Tv_O2 += exp(-(veO2[i]-i*veO2[0])/(KB*T)-
			                i*veO2[0]/(KB*Tv[2]));
    }
    std::cout << "Non-eq. partition functions" << std::endl;
    std::cout << z_vibr_T_Tv_NO << " " 
	      << z_vibr_T_Tv_N2 << " " 	
	      << z_vibr_T_Tv_O2 << std::endl;

    // Compute state-to-state diss/rec coefficients, k_c,diss^d(0)
    //double kcd_n2[lN2], kcr_n2[lN2];
    //double kcd_o2[lO2], kcr_o2[lO2];
    //double kcd_no[lNO], kcr_no[lNO];

    // Mass and diameter
    double m_NO = 4.9826300488143997e-26; 	// kg
    double d_NO = 3.4061e-10; 			// m
    double m_N2 = 4.6517344343135997e-26;
    double d_N2 = 3.40385035355259e-10; 
    double m_O2 = 5.3135256633152e-26;
    double d_O2 = 3.5155e-10; 
    double m_N  = 2.3258672171567998e-26;
    double m_O  = 2.6567628316576e-26;

    double m[5], d[5];
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
    
    double kcd_no[5][lNO]; 
    double coll_mass = 0.;
    double diameter = 0.;
    for (int i=0; i<5; ++i) { // N O NO N2 O2
	coll_mass = (m[i] * m[2]) / (m[i] + m[2]);
	diameter = 0.5 * (d[i] + d[2]);
        for (int j=0; j<lNO; ++j) 
            kcd_no[i][j] = k_diss_RS(T, coll_mass, diameter, deNO, 
	        veNO[j], /*center_of_mass=*/true);
	//std::cout << i << " " << j << " " << kcd_no[i][j] << std::endl;
    }

    double kcd_n2[5][lN2]; 
    for (int i=0; i<5; ++i) { // N O NO N2 O2
	coll_mass = (m[i] * m[3]) / (m[i] + m[3]);
	diameter = 0.5 * (d[i] + d[3]);
        for (int j=0; j<lN2; ++j) 
            kcd_n2[i][j] = k_diss_RS(T, coll_mass, diameter, deN2, 
	        veN2[j], /*center_of_mass=*/true);
	//std::cout << i << " " << j << " " << kcd_n2[i][j] << std::endl;
    }

    double kcd_o2[5][lO2]; 
    for (int i=0; i<5; ++i) { // N O NO N2 O2
	coll_mass = (m[i] * m[4]) / (m[i] + m[4]);
	diameter = 0.5 * (d[i] + d[4]);
        for (int j=0; j<lO2; ++j) 
            kcd_o2[i][j] = k_diss_RS(T, coll_mass, diameter, deO2, 
	        veO2[j], /*center_of_mass=*/true);
	//std::cout << i << " " << j << " " << kcd_o2[i][j] << std::endl;
    }

//    // Compute dissociation coefficients, k_c,diss^d(0)
//    // according to the model of Rigid Sphere (RS).
//    double coll_mass = 0.5 * m_NO; // (m1 * m2)/(m1 + m2)
//    double diameter = d_NO; // 0.5 * (d1 + d2)
//
//    // TODO: instead of calling lNO times the function, pass directly
//    // the whole array of vibrational energy.
//    std::cout << "kcd_no" << std::endl;
//    for (int i=0; i<lNO; ++i) {
//
//        kcd_no[i] = k_diss_RS(T, coll_mass, diameter, deNO, 
//	    veNO[i], /*center_of_mass*/true);
//	std::cout << i << " " << kcd_no[i] << std::endl;
//    }
//
//    coll_mass = 0.5 * m_N2;
//    diameter = d_N2; 
//    std::cout << "kcd_n2" << std::endl;
//    for (int i=0; i<lN2; ++i) {
//
//        kcd_n2[i] = k_diss_RS(T, coll_mass, diameter, deN2, 
//	    veN2[i], /*center_of_mass*/true);
//	std::cout << i << " " << kcd_n2[i] << std::endl;
//    }
//
//    coll_mass = 0.5 * m_O2;
//    diameter = d_O2; 
//    std::cout << "kcd_o2" << std::endl;
//    for (int i=0; i<lO2; ++i) {
//
//        kcd_o2[i] = k_diss_RS(T, coll_mass, diameter, deO2, 
//	    veO2[i], /*center_of_mass*/true);
//	std::cout << i << " " << kcd_o2[i] << std::endl;
//    }

    // Compute the rotational partition function
    // according to a simplified formula ... TODO: better formula
    double Be[3] = {1.998, 1.4377, 1.6720}; // *100 // m^-1
    int sigma[3] = {2, 2, 1};
    double theta_r[3] = {0, 0, 0};
    double z_rot[3] = {0, 0, 0};
    for (int i=0; i<3; ++i) {
        theta_r[i] = Be[i]*100*HP*C0/KB;
        z_rot[i] = T/(sigma[i]*theta_r[i]);
	//std::cout << theta_r[i] << " " << z_rot[i] << std::endl;
    }

    // Compute STS recombination coefficients, k_rec,c^d(0)
    double fac0 = HP*HP*HP * pow(2*PI*KB*T,-1.5);
    double fac1 = (m_NO/(m_N*m_O)) * fac0 * z_rot[0];
    //std::cout << "kcr_no" << std::endl;

    double kcr_no[5][lNO]; 
    for (int i=0; i<5; ++i) {
        for (int j=0; j<lNO; ++j) 
            kcr_no[i][j] = kcd_no[i][j] * exp(-(veNO[j]-deNO)/(KB*T)) * fac1;
	//std::cout << i << " " << kcr_no[i][j] << std::endl;
    }

    fac1 = (m_N2/(m_N*m_N)) * fac0 * z_rot[1];
    //std::cout << "kcr_n2" << std::endl;

    double kcr_n2[5][lN2]; 
    for (int i=0; i<5; ++i) {
        for (int j=0; j<lN2; ++j) 
            kcr_n2[i][j] = kcd_n2[i][j] * exp(-(veN2[j]-deN2)/(KB*T)) * fac1;
	//std::cout << i << " " << kcr_n2[i][j] << std::endl;
    }

    fac1 = (m_O2/(m_O*m_O)) * fac0 * z_rot[2];
    //std::cout << "kcr_o2" << std::endl;
    
    double kcr_o2[5][lO2]; 
    for (int i=0; i<5; ++i) {
        for (int j=0; j<lO2; ++j) 
            kcr_o2[i][j] = kcd_o2[i][j] * exp(-(veO2[j]-deO2)/(KB*T)) * fac1;
	//std::cout << i << " " << kcd_o2[i][j] << std::endl;
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
            fac0 += exp(-((veNO[j]-i*veNO[0])/(KB*T)) -
                                   i*veNO[0]/(KB*Tv[0])) * kcd_no[i][j];
	    fac1 += kcr_no[i][j];
        }
        kd_no[i] = fac0 / z_vibr_T_Tv_NO;
        kr_no[i] = fac1;
    }

    // N2
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
	fac0 = 0., fac1 = 0.; 
        for (int j=0; j<lN2; ++j) {
            fac0 += exp(-((veN2[j]-i*veN2[0])/(KB*T)) -
                                   i*veN2[0]/(KB*Tv[1])) * kcd_n2[i][j];
	    fac1 += kcr_n2[i][j];
        }
        kd_n2[i] = fac0 / z_vibr_T_Tv_N2;
        kr_n2[i] = fac1;
    }

    // O2
    for (int i=0; i<m_thermo.nSpecies(); ++i) {
	fac0 = 0., fac1 = 0.; 
        for (int j=0; j<lO2; ++j) {
            fac0 += exp(-((veO2[j]-i*veO2[0])/(KB*T)) -
                                   i*veO2[0]/(KB*Tv[2])) * kcd_o2[i][j];
	    fac1 += kcr_o2[i][j];
        }
        kd_o2[i] = fac0 / z_vibr_T_Tv_O2;
        kr_o2[i] = fac1;
    }

    //for (int i=0; i<m_thermo.nSpecies(); ++i) {
    //    fac0 = 0., fac1 = 0.;    
    //    for (int j=0; j<lNO; ++j) {
    //        fac0 += exp(-((veNO[i][j]-i*veNO[0])/(KB*T)) - 
    //    	   	              i*veNO[0]/(KB*Tv[0])) * kcd_no[i][j];
    //        fac1 += kcr_no[i][j];
    //    }
    // 	    kd[i][j] = fac0 * z_vibr_T_Tv_NO ; // TODO: reshape z_vibr_T_Tv_NO
    //        kr[i][j] += kcr_no[i][j];
    //}
    //double kNO = tmp * z_vibr_T_Tv_NO;
    //kd[2] = kNO;

    //tmp = 0.;
    //for (int i=0; i<lN2; ++i) {
    //    tmp += exp(-((veN2[i]-i*veN2[0])/(KB*T)) - 
    //		              i*veN2[0]/(KB*Tv[1])) * kcd_n2[i];
    //    kr[3] += kcr_n2[i];
    //}
    //double kN2 = tmp * z_vibr_T_Tv_N2 ;
    //kd[3] = kN2;

    //tmp = 0.;
    //for (int i=0; i<lO2; ++i) {
    //    tmp += exp(-((veO2[i]-i*veO2[0])/(KB*T)) - 
    //    		      i*veO2[0]/(KB*Tv[1])) * kcd_o2[i];
    //    kr[4] += kcr_o2[i];
    //}
    //double kO2 = tmp * z_vibr_T_Tv_O2 ;
    //kd[4] = kO2;

    //std::cout << "kd and kr" << std::endl;
    //std::cout << kd[0] << "   " << kr[0] << "\n"
    //	      << kd[1] << "   " << kr[1] << "\n" 
    //	      << kd[2] << "   " << kr[2] << "\n"  
    //	      << kd[3] << "   " << kr[3] << "\n" 
    //	      << kd[4] << "   " << kr[4] << std::endl;

    // Here, as elsewhere, a certain fixed order of species is assumed.
    // Such order is given by the list of species in the input file:
    // N O NO N2 O2
    double R_react_23[m_thermo.nSpecies()] = {0.,0.,0.,0.,0.};
    double nd[m_thermo.nSpecies()] = {nd_n, nd_o, nd_no, nd_n2, nd_o2};

    // N
    R_react_23[0] = 0.;
    // O
    R_react_23[1] = 0.;
    // NO
    //R_react_23[2] = nd[0] * (nd[0]*nd[1]*kr[0]-nd[2]*kd[0]) +
    //                nd[1] * (nd[0]*nd[1]*kr[1]-nd[2]*kd[1]) +
    //                nd[2] * (nd[0]*nd[1]*kr[2]-nd[2]*kd[2]) +
    //                nd[3] * (nd[0]*nd[1]*kr[3]-nd[2]*kd[3]) +
    //                nd[4] * (nd[0]*nd[1]*kr[4]-nd[2]*kd[4]) ;
    
    R_react_23[2] = nd[0] * (nd[0]*nd[1]*kr_no[0]-nd[2]*kd_no[0]) +
                    nd[1] * (nd[0]*nd[1]*kr_no[1]-nd[2]*kd_no[1]) +
                    nd[2] * (nd[0]*nd[1]*kr_no[2]-nd[2]*kd_no[2]) +
                    nd[3] * (nd[0]*nd[1]*kr_no[3]-nd[2]*kd_no[3]) +
                    nd[4] * (nd[0]*nd[1]*kr_no[4]-nd[2]*kd_no[4]) ;
    // N2
    //r_react_23[3] = nd[0] * (nd[0]*nd[0]*kr[0]-nd[3]*kd[0]) +
    //                nd[1] * (nd[0]*nd[0]*kr[1]-nd[3]*kd[1]) +
    //                nd[2] * (nd[0]*nd[0]*kr[2]-nd[3]*kd[2]) +
    //                nd[3] * (nd[0]*nd[0]*kr[3]-nd[3]*kd[3]) +
    //                nd[4] * (nd[0]*nd[0]*kr[4]-nd[3]*kd[4]) ;
    
    R_react_23[3] = nd[0] * (nd[0]*nd[0]*kr_n2[0]-nd[3]*kd_n2[0]) +
                    nd[1] * (nd[0]*nd[0]*kr_n2[1]-nd[3]*kd_n2[1]) +
                    nd[2] * (nd[0]*nd[0]*kr_n2[2]-nd[3]*kd_n2[2]) +
                    nd[3] * (nd[0]*nd[0]*kr_n2[3]-nd[3]*kd_n2[3]) +
                    nd[4] * (nd[0]*nd[0]*kr_n2[4]-nd[3]*kd_n2[4]) ;
    // O2
    //R_react_23[4] = nd[0] * (nd[1]*nd[1]*kr[0]-nd[4]*kd[0]) +
    //                nd[1] * (nd[1]*nd[1]*kr[1]-nd[4]*kd[1]) +
    //                nd[2] * (nd[1]*nd[1]*kr[2]-nd[4]*kd[2]) +
    //                nd[3] * (nd[1]*nd[1]*kr[3]-nd[4]*kd[3]) +
    //                nd[4] * (nd[1]*nd[1]*kr[4]-nd[4]*kd[4]) ;

    R_react_23[4] = nd[0] * (nd[1]*nd[1]*kr_o2[0]-nd[4]*kd_o2[0]) +
                    nd[1] * (nd[1]*nd[1]*kr_o2[1]-nd[4]*kd_o2[1]) +
                    nd[2] * (nd[1]*nd[1]*kr_o2[2]-nd[4]*kd_o2[2]) +
                    nd[3] * (nd[1]*nd[1]*kr_o2[3]-nd[4]*kd_o2[3]) +
                    nd[4] * (nd[1]*nd[1]*kr_o2[4]-nd[4]*kd_o2[4]) ;

    // TODO: check dimensions!
    for (int i = 0; i < m_thermo.nSpecies(); ++i) {
        p_r_23[i] = R_react_23[i];
	std::cout << p_r_23[i] << std::endl;
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

        double conversion = sqrt(2 * KB * T / coll_mass);

        auto integrand = [T, degree, coll_mass, diameter, diss_energy, 
	    vibr_energy, center_of_mass, conversion](double g)
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
        auto integrand = [T, degree, coll_mass, vss_c_cs, vss_omega, 
            diss_energy, vibr_energy, center_of_mass, conversion](double g) 
	    {
                return pow(g, 2 * degree + 3) * crosssection_diss_VSS(
	            conversion * g, coll_mass, vss_c_cs, vss_omega, 
		        diss_energy, vibr_energy, center_of_mass) * 
			    exp(-g * g); 
	    };

        return sqrt(KB * T / (2 * PI * coll_mass)) * 
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

        double energy = vibr_energy + rel_vel * rel_vel * coll_mass / 2;

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
        return sqrt(2 * (diss_energy - vibr_energy) / coll_mass);
    }

//==============================================================================

    double Kinetics::crosssection_elastic_RS(double diameter) 
    {
        return PI * diameter * diameter;
    }

    double Kinetics::crosssection_elastic_VSS(double rel_vel, double vss_c_cs, 
	double vss_omega) 
    {
        return vss_c_cs * pow(rel_vel, 1. - 2. * vss_omega);
    }

    double Kinetics::crosssection_elastic_VSS(double rel_vel, 
	double coll_mass, double vss_c, double vss_omega) 
    {
        return vss_c * pow(coll_mass * rel_vel * rel_vel / (2 * KB), 
	    -vss_omega);
    }

//==============================================================================

    double Kinetics::k_Arrhenius(double T, double arrhenius_A, 
	double arrhenius_n, double energy) 
    {
        return arrhenius_A * pow(T, arrhenius_n) * exp(-energy / (KB * T));
    }

//==============================================================================

//   double Kinetics::Z_diss(double T, const double * vib_energy, 
//       int num_vibr_levels, int i) 
//   {
//       return p_Z_vibr_eq(T, vibr_energy) * 
//	   exp(vibr_energy[i] / (KB * T)) / num_vibr_levels;
//   }

//==============================================================================
  
//   double Kinetics::p_Z_vibr_eq(double T, const double *vibr_energy) 
//   {
//       double p_Z_vibr_eq = 0.;
//       for (int i=0; i<vibr_energy.size(); ++i)
//           p_Z_vibr_eq += exp(-vibr_energy[i] / (KB * T));	       
//       return p_Z_vibr_eq;
//   }

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
