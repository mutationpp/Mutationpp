/**
 * @file ReactionType.h
 *
 * @brief Defines ReationType enumeration.
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

#ifndef KINETICS_REACTION_TYPE_H
#define KINETICS_REACTION_TYPE_H


namespace Mutation {
    namespace Kinetics {

/**
 * Enumerates possible reaction types.
 * @see Reaction::type()
 */    
enum ReactionType
{
    ASSOCIATIVE_IONIZATION,
    CHARGE_EXCHANGE,
    DISSOCIATION_E,
    DISSOCIATION_M,
    DISSOCIATIVE_RECOMBINATION,
    ELECTRONIC_ATTACHMENT,
    ELECTRONIC_DETACHMENT,
    EXCHANGE,
    IONIZATION_E,
    IONIZATION_M,
    ION_RECOMBINATION_E,
    ION_RECOMBINATION_M,
    RECOMBINATION_E,
    RECOMBINATION_M,
    EXCITATION_E,
    EXCITATION_M,
    DEEXCITATION_E,
    DEEXCITATION_M,
    MAX_REACTION_TYPES /// Number of ReactionType values (leave at end of list)
};

/**
 * Simpy returns a string which represents the reaction type given by the 
 * ReactionType argument.
 */
const char* const reactionTypeString(const ReactionType type);

    } // namespace Kinetics
} // namespace Mutation

#endif // KINETICS_REACTION_TYPE_H
