#ifndef SETUPLARSENTTv_H
#define SETUPLARSENTTv_H

#include "SetupProblem.h"

class SetupLarsenTTv: public SetupProblem {
public:
    SetupLarsenTTv();
    ~SetupLarsenTTv(){}

    Data* getData(Mutation::Mixture& l_mix);
    Problem* getProblem(Mutation::Mixture& l_mix, Data& l_data);
};

#endif /* SETUPLARSENTTv_H */
