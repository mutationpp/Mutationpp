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
