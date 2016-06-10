#ifndef RANKINEHUGONIOTNT_H
#define RANKINEHUGONIOTNT_H

#include "mutation++.h"

#include "ShockRelations.h"

class RankineHugoniotNT: public ShockRelations {
public:
    RankineHugoniotNT(Mutation::Mixture& l_mix);
    ~RankineHugoniotNT();

    void applyShockRelations(const Data& l_data_pre, Data& l_data_post);

private:
    Mutation::Mixture& m_mix;
    const size_t set_state_rhoi_T;
};

#endif /* RANKINEHUGONIOT1T_H */
