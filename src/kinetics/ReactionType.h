/**
 * @file ReactionType.h
 *
 * @brief Defines ReationType enumeration.
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
    DISSOCIATIVE_RECOMBINATION,
    ASSOCIATIVE_DETACHMENT,
    DISSOCIATIVE_ATTACHMENT,
    DISSOCIATION_E,
    RECOMBINATION_E,
    DISSOCIATION_M,
    RECOMBINATION_M,
    IONIZATION_E,
    ION_RECOMBINATION_E,
    IONIZATION_M,
    ION_RECOMBINATION_M,
    ELECTRONIC_DETACHMENT_E,
    ELECTRONIC_ATTACHMENT_E,
    ELECTRONIC_DETACHMENT_M,
    ELECTRONIC_ATTACHMENT_M,
    EXCHANGE,
    EXCITATION_E,
    EXCITATION_M,
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
