/**
 * @file ThermoDB.cpp
 *
 * @brief Implementation of the ThermoDB class.
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

#include "ThermoDB.h"
#include "Utilities.h"

#include <Eigen/Dense>
using namespace Eigen;

#include <algorithm>
#include <cassert>
using namespace std;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ThermoDB::ThermoDB(double sst, double ssp)
    : m_sst(sst), m_ssp(ssp)
{
    assert(m_sst >= 0.0);
    assert(m_ssp >  0.0);
}

//==============================================================================

bool ThermoDB::load(const SpeciesListDescriptor& descriptor)
{
    // It is possible that this isn't the first call to load so just make sure
    // species and elements are cleared
    m_species.clear();
    m_elements.clear();
    
    // Load all possible species from the concrete database type
    std::list<Species> species_list;
    loadAvailableSpecies(species_list);
    
    // Check for duplicate species in the database
    std::list<Species>::iterator iter1 = species_list.begin();
    std::list<Species>::iterator iter2;
    
    while (iter1 != species_list.end()) {
        (iter2 = iter1)++;
        while (iter2 != species_list.end()) {
            if (iter1->name() == iter2->name()) {
                std::cout << "Warning, species \"" << iter1->name()
                          << "\" is defined more than once in the thermodynamic"
                          << " database.  I will ignore recurrences..."
                          << std::endl;
                iter2 = species_list.erase(iter2);
                while (iter2 != species_list.end()) {
                    if (iter1->name() == iter2->name())
                        iter2 = species_list.erase(iter2);
                    else
                        iter2++;
                }
            } else
                iter2++;
        }
        iter1++;
    }
    
    // Use the species list descriptor to remove unwanted species
    iter1 = species_list.begin();
    while (iter1 != species_list.end()) {
        if (descriptor.matches(*iter1))
            iter1++;
        else
            iter1 = species_list.erase(iter1);
    }
    
    // Now we have all of the species that we want but possibly in the wrong
    // order so use the descriptor to tell us the correct order
    std::vector<std::string> missing;
    descriptor.order(species_list, m_species, missing);
    
    if (missing.size() > 0) {
        std::cout << "Missing the following species!" << std::endl;
        for (int i = 0; i < missing.size(); ++i)
            std::cout << "  " << missing[i] << std::endl;
        return false;
    }
    
    // Finally fill our elements vector with only the elements that are required
    // for the species list
    std::set<std::string> element_names;
    for (int i = 0; i < m_species.size(); ++i) {
        Species::StoichList::const_iterator iter =
            m_species[i].stoichiometry().begin();
        for ( ; iter != m_species[i].stoichiometry().end(); ++iter)
            element_names.insert(iter->first);
    }
    
    for (int i = 0; i < Element::database().size(); ++i)
        if (element_names.count(Element::database()[i].name()) > 0)
            m_elements.push_back(Element::database()[i]);
    
    // Tell the concrete class to load the necessary thermodynamic data based on
    // the final species list
    loadThermodynamicData();
    return true;
}

//==============================================================================

// Used in ThermoDB::cv
struct MinusOne {
    void operator () (double& x) const { x -= 1.0; }
} MinusOne;

void ThermoDB::cv(
    double Th, double Te, double Tr, double Tv, double Tel, 
    double* const cv = NULL, double* const cvt = NULL, double* const cvr = NULL,
    double* const cvv = NULL, double* const cvel = NULL)
{
    // Compute Cp/Ru
    cp(Th, Te, Tr, Tv, Tel, cv, cvt, cvr, cvv, cvel);
    
    const size_t ns = m_species.size();

    // Cv/Ru = Cp/Ru - 1
    if (cv   != NULL) std::for_each(  cv+0,   cv+ns, MinusOne);
    if (cvt  != NULL) std::for_each( cvt+0,  cvt+ns, MinusOne);
    if (cvr  != NULL) std::for_each( cvr+0,  cvr+ns, MinusOne);
    if (cvv  != NULL) std::for_each( cvv+0,  cvv+ns, MinusOne);
    if (cvel != NULL) std::for_each(cvel+0, cvel+ns, MinusOne);
}

//==============================================================================

void ThermoDB::cpv(double T, double* const p_cp)
{
    throw NotImplementedError("ThermoDB::cpv()");
}

//==============================================================================

void ThermoDB::cpel(double T, double* const p_cp)
{
    throw NotImplementedError("ThermoDB::cpel()");
}

//==============================================================================

void ThermoDB::hv(double T, double* const p_h)
{
    throw NotImplementedError("ThermoDB::hv()");
}

//==============================================================================

void ThermoDB::hel(double T, double* const p_h)
{
    throw NotImplementedError("ThermoDB::hel()");
}

//==============================================================================

void ThermoDB::cpint(double T, double* const p_cp) {
    cp(T, T, T, T, T, p_cp);
    Map<ArrayXd>(p_cp, m_species.size()) -= 2.5;
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation

