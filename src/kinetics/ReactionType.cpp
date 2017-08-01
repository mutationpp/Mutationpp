/**
 * @file ReactionType.cpp
 *
 * @brief Implementation of ReationType related functions.
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

#include "ReactionType.h"

namespace Mutation {
    namespace Kinetics {

//==============================================================================

const char* const reactionTypeString(const ReactionType type)
{
    switch (type) {
        case ASSOCIATIVE_IONIZATION:
            return "associative ionization";
        case CHARGE_EXCHANGE:
            return "charge exchange";
        case DISSOCIATION_M:
            return "heavy particle impact dissociation";
        case DISSOCIATION_E:
            return "electron impact dissociation";
        case RECOMBINATION_M:
            return "heavy particle impact recombination";
        case RECOMBINATION_E:
            return "electron impact recombination";
        case DISSOCIATIVE_RECOMBINATION:
            return "dissociative recombination";
        case ELECTRONIC_ATTACHMENT:
            return "electron attachment";
        case ELECTRONIC_DETACHMENT:
            return "electron detachment";
        case EXCHANGE:
            return "exchange";
        case IONIZATION_M:
            return "heavy particle impact ionization";
        case IONIZATION_E:
            return "electron impact ionization";
        case ION_RECOMBINATION_M:
            return "heavy particle impact ion recombination";
        case ION_RECOMBINATION_E:
            return "electron impact ion recombination";
	case EXCITATION_E:
 	    return "electron impact excitation";
	case EXCITATION_M:
 	    return "heavy particle impact excitation";	
	case DEEXCITATION_E:
 	    return "electron impact deexcitation";
	case DEEXCITATION_M:
 	    return "heavy particle impact deexcitation";
        default:
            return "unknown type";
    }
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation
