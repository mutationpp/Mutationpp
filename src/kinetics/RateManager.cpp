#include <iostream>
#include <typeinfo>

#include "RateManager.h"
#include "Reaction.h"
#include "StateModel.h"

namespace Mutation {
    namespace Kinetics {

//==============================================================================
    
// Simple macro to create a temperature selector type
#define TEMPERATURE_SELECTOR(__NAME__,__T__)\
class __NAME__\
{\
public:\
    inline double getT(const Thermodynamics::StateModel* const state) const {\
        return ( __T__ );\
    }\
};

/// Temperature selector which returns the current translational temperature
TEMPERATURE_SELECTOR(TSelector, state->T())

/// Temperature selector which returns the current electron temperature
TEMPERATURE_SELECTOR(TeSelector, state->Te())

/// Temperature selector which returns the current value of sqrt(T*Tv)
TEMPERATURE_SELECTOR(ParkSelector, std::sqrt(state->T()*state->Tv()))

#undef TEMPERATURE_SELECTOR

/// Arrhenius group evaluated at T
typedef RateLawGroup1T<Arrhenius, TSelector> ArrheniusT;

/// Arrhenius group evaluated at Te
typedef RateLawGroup1T<Arrhenius, TeSelector> ArrheniusTe;

/// Arrhenius group evaluated at sqrt(T*Tv)
typedef RateLawGroup1T<Arrhenius, ParkSelector> ArrheniusPark; 

//==============================================================================

/**
 * Used to define which rate law groups the forward and reverse rate laws are
 * evaluated in for a given ReactionType value.  The default is an Arrhenius
 * rate law with Tf = Tb = T.
 */
template <int Type>
struct RateSelector {
    typedef ArrheniusT ForwardGroup;
    typedef ArrheniusT ReverseGroup;
};

#define SELECT_RATE_LAWS(__TYPE__,__FORWARD__,__REVERSE__)\
template <> struct RateSelector<__TYPE__> {\
    typedef __FORWARD__ ForwardGroup;\
    typedef __REVERSE__ ReverseGroup;\
};

// Default rate law groups for non (kf(T), kb(T)) reaction types
SELECT_RATE_LAWS(ASSOCIATIVE_IONIZATION, ArrheniusT, ArrheniusTe)
SELECT_RATE_LAWS(DISSOCIATION_E, ArrheniusTe, ArrheniusTe)
SELECT_RATE_LAWS(DISSOCIATION_M, ArrheniusPark, ArrheniusT)
SELECT_RATE_LAWS(DISSOCIATIVE_RECOMBINATION, ArrheniusTe, ArrheniusT)
SELECT_RATE_LAWS(ELECTRONIC_ATTACHMENT, ArrheniusTe, ArrheniusT)
SELECT_RATE_LAWS(ELECTRONIC_DETACHMENT, ArrheniusT, ArrheniusTe)
SELECT_RATE_LAWS(IONIZATION_E, ArrheniusTe, ArrheniusTe)
SELECT_RATE_LAWS(IONIZATION_M, ArrheniusT, ArrheniusTe)
SELECT_RATE_LAWS(ION_RECOMBINATION_E, ArrheniusTe, ArrheniusTe)
SELECT_RATE_LAWS(ION_RECOMBINATION_M, ArrheniusTe, ArrheniusT)
SELECT_RATE_LAWS(RECOMBINATION_E, ArrheniusTe, ArrheniusTe)
SELECT_RATE_LAWS(RECOMBINATION_M, ArrheniusT, ArrheniusPark)

#undef SELECT_RATE_LAWS

//==============================================================================

RateManager::RateManager(const std::vector<Reaction>& reactions)
    : m_nr(reactions.size())
{
    // Add all of the reactions' rate coefficients to the manager
    const size_t nr = reactions.size();
    for (size_t i = 0; i < nr; ++i)
        addReaction(i, reactions[i]);
    
    // Allocate storage in one block for both rate coefficient arrays
    mp_lnkff = new double [2*nr];
    mp_lnkfb = mp_lnkff + nr;
    
    // Initialize the arrays to zero
    std::fill(mp_lnkff, mp_lnkff+2*nr, 0.0);
}

//==============================================================================

RateManager::~RateManager()
{
    // Storage for lnkff and lnkfb were both allocated in mp_lnkff
    if (mp_lnkff != NULL)
        delete [] mp_lnkff;
}

//==============================================================================

void RateManager::addReaction(const size_t rxn, const Reaction& reaction)
{
    // Get the rate law which is being used in this reaction
    const RateLaw* p_rate = reaction.rateLaw();
    
    // Arrhenius reactions
    if (typeid(*p_rate) == typeid(Arrhenius)) {
        selectRate<MAX_REACTION_TYPES-1>(reaction.type(), rxn, p_rate);
    } else {
        std::cerr << "Rate law " << typeid(*p_rate).name()
                  << " not implemented in RateManager!" << std::endl;
        exit(1);
    }
    
    cout << "Number of reactions with Tf = Tb: " << m_to_copy.size() << endl;
}

//==============================================================================

template <int NReactionTypes>
void RateManager::selectRate(
    const ReactionType type, const size_t rxn, const RateLaw* const p_rate)
{
    if (type == NReactionTypes)
        addRate<
            typename RateSelector<NReactionTypes>::ForwardGroup,
            typename RateSelector<NReactionTypes>::ReverseGroup>(rxn, p_rate);
    else
        selectRate<NReactionTypes-1>(type, rxn, p_rate);
}

template <>
void RateManager::selectRate<0>(
    const ReactionType type, const size_t rxn, const RateLaw* const p_rate)
{
    addRate<
        typename RateSelector<0>::ForwardGroup,
        typename RateSelector<0>::ReverseGroup>(rxn, p_rate);
}

//==============================================================================

template <typename T, typename U>
struct is_same {
    enum { value = 0 };
};

template <typename T>
struct is_same<T,T> {
    enum { value = 1 };
};

//==============================================================================

template <typename ForwardGroup, typename ReverseGroup>
void RateManager::addRate(const size_t rxn, const RateLaw* const p_rate)
{    
    m_rate_groups.addRateCoefficient<ForwardGroup>(rxn, p_rate);
    
    /// Make use of forward computation when possible
    if (is_same<ForwardGroup, ReverseGroup>::value)
        m_to_copy.push_back(rxn);
    else
        // Evaluate at the reverse temperature
        // note: mp_lnkff+(rxn+m_nr) = mp_lnkfb+rxn
        m_rate_groups.addRateCoefficient<ReverseGroup>(rxn+m_nr, p_rate);
}

//==============================================================================

void RateManager::update(const Thermodynamics::StateModel* p_state)
{
    m_rate_groups.logOfRateCoefficients(p_state, mp_lnkff);
    
    std::vector<size_t>::const_iterator iter = m_to_copy.begin();
    for ( ; iter != m_to_copy.end(); ++iter)
        mp_lnkfb[*iter] = mp_lnkff[*iter];
}

//==============================================================================

    } // namespace Kinetics
} // namespace Mutation

