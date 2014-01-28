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
            return "atom impact excitation";
        case DEEXCITATION_E:
            return "electron impact de-excitation";
        case DEEXCITATION_M:
            return "atom impact de-excitation";
        default:
            return "unknown type";
    }
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation
