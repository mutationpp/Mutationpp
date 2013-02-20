#include <iostream>
#include <cstdlib>

#include "StoichiometryManager.h"

using namespace std;
using namespace Mutation::Numerics;

namespace Mutation {
    namespace Kinetics {

void StoichiometryManager::addReaction(const int rxn, const vector<size_t> &sps)
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
            std::cout << "Error trying to add reaction with more than 3\n";
            std::cout << "species on a single side!" << std::endl;
            exit(1);
    }
}

#define STOICH_MGR_APPLY_FUNC(__my_func__,__stoic_func__,__in__,__out__)\
template<class Iterator, class Vec1, class Vec2>\
inline static void _##__my_func__ (\
    Iterator begin, const Iterator end, const Vec1& input, Vec2& output)\
{\
    for (; begin != end; ++begin)\
        begin-> __stoic_func__ (input, output);\
}\
void StoichiometryManager:: __my_func__ (\
    const RealVector& __in__ , RealVector& __out__ ) const\
{\
    _##__my_func__ (m_stoich1_vec.begin(), m_stoich1_vec.end(), __in__ , __out__ );\
    _##__my_func__ (m_stoich2_vec.begin(), m_stoich2_vec.end(), __in__ , __out__ );\
    _##__my_func__ (m_stoich3_vec.begin(), m_stoich3_vec.end(), __in__ , __out__ );\
}

STOICH_MGR_APPLY_FUNC(multReactions, multReaction, s, r)
STOICH_MGR_APPLY_FUNC(incrReactions, incrReaction, s, r)
STOICH_MGR_APPLY_FUNC(decrReactions, decrReaction, s, r)

STOICH_MGR_APPLY_FUNC(incrSpecies, incrSpecies, r, s)
STOICH_MGR_APPLY_FUNC(decrSpecies, decrSpecies, r, s)

#undef STOICH_MGR_APPLY_FUNC


    } // namespace Kinetics
} // namespace Mutation

