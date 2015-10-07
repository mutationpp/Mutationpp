#include <iostream>
#include <cstdlib>

#include "GSIStoichiometryManager.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//========================================================================================

void GSIStoichiometryManagerGamma::addReaction(const int rxn, const std::vector<int>& sps){
    switch (sps.size()){
        case 1:
            m_stoich1_vec.push_back(Stoich1(rxn, sps[0]));
            break;
        case 2:
            m_stoich2_vec.push_back(Stoich2(rxn, sps[0], sps[1]));
            break;
        default:
            std::cout << "Error trying to add Gas Surface interaction reaction with\n";
            std::cout << "more than 2 species on a single side with the gamma model!" << std::endl;
            exit(1);
    }
}

//========================================================================================

#define STOICH_MGR_APPLY_FUNC(__my_func__,__stoic_func__)\
template<class Iterator, class Vec1, class Vec2>\
inline static void _##__my_func__ (\
    Iterator begin, const Iterator end, const Vec1& input, Vec2& output)\
{\
    for (; begin != end; ++begin)\
        begin-> __stoic_func__ (input, output);\
}\
void GSIStoichiometryManagerGamma:: __my_func__ (\
    const Eigen::VectorXd& in, Eigen::VectorXd& out) const\
{\
    _##__my_func__ (m_stoich1_vec.begin(), m_stoich1_vec.end(), in , out );\
    _##__my_func__ (m_stoich2_vec.begin(), m_stoich2_vec.end(), in , out );\
}\

STOICH_MGR_APPLY_FUNC(multReactions, multReaction)
STOICH_MGR_APPLY_FUNC(incrReactions, incrReaction)
STOICH_MGR_APPLY_FUNC(decrReactions, decrReaction)

STOICH_MGR_APPLY_FUNC(incrSpecies, incrSpecies)
STOICH_MGR_APPLY_FUNC(decrSpecies, decrSpecies)

#undef STOICH_MGR_APPLY_FUNC

//=======================================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
