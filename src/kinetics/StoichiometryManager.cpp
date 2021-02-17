/**
 * @file StoichiometryManager.cpp
 *
 * @brief Implementation of StoichiometryManager class.
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

#include <iostream>
#include <cstdlib>

#include "Errors.h"
#include "StoichiometryManager.h"

using namespace std;

namespace Mutation {
    namespace Kinetics {

void StoichiometryManager::addReaction(const int rxn, const vector<int>& sps)
{
    switch (sps.size()) {
        case 1:
            m_stoich1_vec.push_back(Stoich1(rxn, sps[0]));
            break;
        case 2:
            m_stoich2_vec.push_back(Stoich2(rxn, sps[0], sps[1]));
            break;
        case 3:
            m_stoich3_vec.push_back(Stoich3(rxn, sps[0], sps[1], sps[2]));
            break;
        default:
            throw InvalidInputError("number of species", sps.size())
                << "Error trying to add reaction with more than 3 "
                << "species on a single side.";
    }
}

#define STOICH_MGR_APPLY_FUNC(__my_func__,__stoic_func__)\
template<class Iterator, class Vec1, class Vec2>\
inline static void _##__my_func__ (\
    Iterator begin, const Iterator end, const Vec1& input, Vec2& output)\
{\
    for (; begin != end; ++begin)\
        begin-> __stoic_func__ (input, output);\
}\
void StoichiometryManager:: __my_func__ (\
    const double* const in, double* const out) const\
{\
    _##__my_func__ (m_stoich1_vec.begin(), m_stoich1_vec.end(), in , out );\
    _##__my_func__ (m_stoich2_vec.begin(), m_stoich2_vec.end(), in , out );\
    _##__my_func__ (m_stoich3_vec.begin(), m_stoich3_vec.end(), in , out );\
}

STOICH_MGR_APPLY_FUNC(multReactions, multReaction)
STOICH_MGR_APPLY_FUNC(incrReactions, incrReaction)
STOICH_MGR_APPLY_FUNC(decrReactions, decrReaction)

STOICH_MGR_APPLY_FUNC(incrSpecies, incrSpecies)
STOICH_MGR_APPLY_FUNC(decrSpecies, decrSpecies)

#undef STOICH_MGR_APPLY_FUNC


    } // namespace Kinetics
} // namespace Mutation

