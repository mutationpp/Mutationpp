#include "RankineHugoniot1T.h"

RankineHugoniot1T::RankineHugoniot1T(Mutation::Mixture& l_mix)
                                        : m_mix(l_mix),
                                          set_state_rhoi_T(1){}

RankineHugoniot1T::~RankineHugoniot1T(){}

void RankineHugoniot1T::applyShockRelations(const Data& l_data_before, Data& l_data_after){

    double l_u1 = l_data_before.getVelocity();
    double l_p1 = l_data_before.getPressure();
    double l_T1 = l_data_before.getTTrans();

    // Set the state with rhoi, T
    m_mix.setState(&l_data_before.getPartialDensities()[0],
                   &l_T1,
                   set_state_rhoi_T);

    double l_gamma =  m_mix.mixtureFrozenCpMass()/m_mix.mixtureFrozenCvMass();
    double l_c1 = m_mix.frozenSoundSpeed();

    double l_gp1 = l_gamma+1.;
    double l_gm1 = l_gamma-1.;
    double l_M1 = l_u1/l_c1;
    double l_M1s = l_M1 * l_M1;

    l_data_after.setPressure(l_p1 * (2.0 * l_gamma * l_M1s - l_gm1) / l_gp1);
    l_data_after.setVelocity(l_u1 - l_c1 * 2.0 / l_gp1 * (l_M1 - 1.0 / l_M1));
    l_data_after.setTTrans(l_T1 * (2.0 * l_gamma * l_M1s - l_gamma + 1) * (l_gm1 + 2.0 / l_M1s) / (l_gp1 * l_gp1) );
    l_data_after.setMassFractions(l_data_before.getMassFractions());

    l_data_after.buildStatePostShock();

}
