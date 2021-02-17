/**
 * @file Reaction.cpp
 *
 * @brief Implementation of Reaction class.
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

#include <cstdlib>
#include <iostream>
#include <cassert>

#include "Reaction.h"

using namespace Mutation::Utilities;
using namespace Mutation::Thermodynamics;

namespace Mutation {
namespace Kinetics {

//==============================================================================

void swap(Reaction& left, Reaction& right) {
    using std::swap;
    swap(left.m_formula,     right.m_formula);
    swap(left.m_reactants,   right.m_reactants);
    swap(left.m_products,    right.m_products);
    swap(left.m_reversible,  right.m_reversible);
    swap(left.m_thirdbody,   right.m_thirdbody);
    swap(left.m_conserves,   right.m_conserves);
    swap(left.m_thirdbodies, right.m_thirdbodies);
    swap(left.m_type,        right.m_type);
    swap(left.mp_rate,       right.mp_rate);
}

//==============================================================================

Reaction::Reaction(const IO::XmlElement& node, const class Thermodynamics& thermo)
    : m_formula(""),
      m_reversible(true),
      m_thirdbody(false),
      m_conserves(true),
      mp_rate(NULL)
{
    // Make sure this is a reaction type XML element
    assert( node.tag() == "reaction" );

    // Store the reaction formula (must have)
    node.getAttribute("formula", m_formula,
        "No formula specied with reaction!");

    // Parse the formula to determine which species are involved, whether or
    // not this is a third-body reaction, and reversibility of the reaction
    parseFormula(node, thermo);

    // Now loop through the children of this node to determine the other
    // attributes of the reaction
    IO::XmlElement::const_iterator iter = node.begin();
    for ( ; iter != node.end(); ++iter) {
        if (iter->tag() == "arrhenius") {
            mp_rate = new Arrhenius(*iter, order());
        } else if (iter->tag() == "M") {
            if (m_thirdbody) {
                std::vector<std::string> tokens;
                String::tokenize(iter->text(), tokens, ":, ");
                for (int i = 0; i < tokens.size(); i+=2) {
                    int index = thermo.speciesIndex(tokens[i]);
                    if (index == 0 and thermo.hasElectrons()) {
                        iter->parseError(
                            "Electron cannot be a thirdbody species!");
                    }

                    if (index >= 0) {
                        m_thirdbodies.push_back(
                            std::make_pair(index, atof(tokens[i+1].c_str())));
                    } else {
                        iter->parseError((
                            std::string("Thirdbody species ") + tokens[i] +
                            std::string(" is not in mixture list!")).c_str());
                    }
                }
            } else {
                iter->parseError((
                    std::string("This reaction is not a thirdbody reaction") +
                    std::string(" but thirdbodies are given!")).c_str());
            }
        }
    }

    // Make sure we got a RateLaw out of all that
    if (mp_rate == NULL)
        node.parseError("A rate law must be supplied with this reaction!");

    // Check for charge and mass conservation
    const size_t ne = thermo.nElements();
    int sums [ne];
    std::fill(sums, sums+ne, 0.0);
    for (int i = 0; i < nReactants(); ++i)
        for (int k = 0; k < ne; ++k)
            sums[k] += thermo.elementMatrix()(m_reactants[i],k);
    for (int i = 0; i < nProducts(); ++i)
        for (int k = 0; k < ne; ++k)
            sums[k] -= thermo.elementMatrix()(m_products[i],k);
    for (int i = 0; i < ne; ++i)
        m_conserves &= (sums[i] == 0);

    // Figure out what type of reaction this is
    determineType(thermo);
}

//==============================================================================

bool Reaction::operator == (const Reaction& r)
{
    // Either both reactions are thirdbody or not
    if (isThirdbody() == r.isThirdbody()) {
        // A + B <=> AB  =  A + B <=> AB
        // A + B  => AB  =  A + B <=> AB
        // A + B <=> AB  =  A + B  => AB
        // A + B  => AB  =  A + B  => AB
        if (m_reactants == r.m_reactants &&
            m_products == r.m_products)
            return true;

        // A + B <=> AB  =  AB <=> A + B,
        // A + B  => AB  =  AB <=> A + B,
        // A + B <=> AB  =  AB  => A + B
        if (m_reactants == r.m_products &&
            m_products == r.m_reactants &&
            !(isReversible() && r.isReversible()))
            return true;
    }

    // In the case where one reaction is thirdbody, we have to check if the non-
    // thirdbody reaction is a special case of the thirdbody reaction
    // (ie: N + O + M <=> NO + M (alpha_N != 0)  and  2N + O <=> NO + N)

    return false;
}

//==============================================================================

void Reaction::parseFormula(
    const IO::XmlElement& node, const class Thermodynamics& thermo)
{
    // First step is to split the formula into reactant and product
    // strings and determine reversibility of the reaction
    size_t pos = m_formula.find("=");
    if (pos == std::string::npos)
        node.parseError((
            std::string("Reaction formula ") + m_formula +
            std::string(" does not have '=' or '=>'!")).c_str());

    std::string reactants = m_formula.substr(0, pos);
    std::string products;

    if (m_formula[pos+1] == '>') {
        m_reversible = false;
        products = m_formula.substr(pos+2, m_formula.length()-pos-1);
    } else {
        m_reversible = true;
        products = m_formula.substr(pos+1, m_formula.length()-pos);
    }

    // Now that we have reactant and product strings, we can parse each
    // separately using the same algorithm
    parseSpecies(m_reactants, reactants, node, thermo);
    parseSpecies(m_products,  products,  node, thermo);
}

//==============================================================================

void Reaction::parseSpecies(
    std::vector<int>& species, std::string& str, const IO::XmlElement& node,
    const class Thermodynamics& thermo)
{
    // State-Machine states for parsing a species formula
    enum ParseState {
        coefficient,
        name,
        plus
    } state = coefficient;

    size_t c = 0;
    size_t s = 0;
    size_t e = 0;
    int nu   = 1;
    bool add_species = false;

    // Remove all spaces to simplify the state machine logic
    String::eraseAll(str, " ");

    // Loop over every character
    while (c != str.size()) {
        // Decide what action to do
        switch(state) {
            case coefficient:
                if (isdigit(str[c])) {
                    nu = atoi(str.substr(c,1).c_str());
                    s = c + 1;
                } else {
                    nu = 1;
                    s = c;
                }
                state = name;
                break;
            case name:
                if (str[c] == '+')
                    state = plus;
                break;
            case plus:
                if (str[c] == '(')  // Allow ionized STS: ie Ar+(0)
                    state = name;
                else if (str[c] != '+') {
                    e = c - 2;
                    c--;
                    add_species = true;
                    state = coefficient;
                }
                break;
        }

        // Add the last species which would not be counted because c will
        // equal str.size() on the next iteration
        if (c == str.size() - 1) {
            add_species = true;
            e = c;
        }

        // If we found the start and end position of a species, add it to
        // the list nu times (unless it is the special case of 'M')
        if (add_species) {
            if (str.substr(s,e-s+1) == "M") {
                m_thirdbody = true;
            } else {
                int index = thermo.speciesIndex(str.substr(s,e-s+1));
                if (index >= 0)
                    for (int i = 0; i < nu; ++i)
                        species.push_back(index);
                else
                    node.parseError(("Species " + str.substr(s,e-s+1) +
                        " is not in the mixture list!").c_str());
            }
            add_species = false;
        }

        // Move on to the next character
        c++;
    }

    // Sort the species vector
    std::sort(species.begin(), species.end());
}

//==============================================================================

void Reaction::determineType(const class Thermodynamics& thermo)
{

    bool reactant_ion      = false;
    bool reactant_electron = false;
    bool product_ion       = false;
    bool product_electron  = false;
    bool m_inert           = isThirdbody();
    bool m_inert_e         = false;
    int ecount             = 0;
    int chgcount           = 0;
    
    // Check reactants for heavy ions and electrons
    for (int i=0; i < nReactants(); ++i) {
        const Species& species = thermo.species(m_reactants[i]);
        if (species.type() == ELECTRON) {
            ecount++;
            reactant_electron = true;
        } else if (species.isIon()) {
            chgcount += abs(species.charge());
            reactant_ion = true;
        }
    }

    // Check products
    for (int i=0; i < nProducts(); ++i) {
        const Species& species = thermo.species(m_products[i]);
        if (species.type() == ELECTRON) {
            ecount--;
            product_electron = true;
        } else if (species.isIon()) {
            chgcount -= abs(species.charge());
            product_ion = true;
        }
    }

    // Check for an inert species if not explicitly a thirdbody reaction
    if (!m_inert) {
        for (int i=0; i < nReactants(); ++i) {
            for (int j=0; j < nProducts(); ++j) {
                m_inert |= (m_reactants[i] == m_products[j]);
            }
        }
    }

    // Is there an inert electron?
    m_inert_e = (product_electron && reactant_electron);

    // ---------------- Logic tree ----------------
    if (ecount > 0) {
        
        if (chgcount > 0) {
           if (m_inert_e)
              m_type = ION_RECOMBINATION_E;
           else if (m_inert)
              m_type = ION_RECOMBINATION_M;
           else 
              m_type = DISSOCIATIVE_RECOMBINATION;
           
        } else {
            if (m_inert_e)
               m_type = ELECTRONIC_ATTACHMENT_E;
            else if (m_inert)
               m_type = ELECTRONIC_ATTACHMENT_M;
            else 
               m_type = DISSOCIATIVE_ATTACHMENT;
        }
        
    } else if (ecount < 0) {
        
        if (chgcount > 0) {
            if (m_inert_e)
                m_type = ELECTRONIC_DETACHMENT_E;
            else if (m_inert)
                m_type = ELECTRONIC_DETACHMENT_M;
            else 
                m_type = ASSOCIATIVE_DETACHMENT;
            
        } else {
            if (m_inert_e)
                m_type = IONIZATION_E;
            else if (m_inert)
                m_type = IONIZATION_M;
            else 
                m_type = ASSOCIATIVE_IONIZATION;
        }
        
    } else {
        
        if (nReactants() > nProducts()) {
            if (m_inert_e)
                m_type = RECOMBINATION_E;
            else
                m_type = RECOMBINATION_M;   // including double recombination        
                
        } else if (nReactants() < nProducts()) {
            if (m_inert_e)
                m_type = DISSOCIATION_E;
            else
                m_type = DISSOCIATION_M;   // including double dissociation  
            
        } else {
            if (m_inert_e)
                m_type = EXCITATION_E;
            else if (m_inert)
                m_type = EXCITATION_M;
            else
                m_type = EXCHANGE;  // exchange of charge, atom, or both simultaneously
        } 
    }
}

//==============================================================================

} // namespace Kinetics
} // namespace Mutation

