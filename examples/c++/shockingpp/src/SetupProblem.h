#ifndef SETUPPROBLEM_H
#define SETUPPROBLEM_H

#include <string>
#include "mutation++.h"

#include "Data.h"
#include "ShockRelations.h"
#include "Problem.h"

class SetupProblem {
public:
    SetupProblem(){}
    virtual ~SetupProblem(){}

    static SetupProblem* createFactory(std::string l_problem_type, std::string l_state_model);

    virtual Data* getDataPreShock(Mutation::Mixture& l_mix, const std::string& l_free_stream_conditions) = 0;
    virtual Data* getDataPostShock(Mutation::Mixture& l_mix) = 0;
    virtual ShockRelations* getShockRelations(Mutation::Mixture& l_mix) = 0;
    virtual Problem* getProblem(Mutation::Mixture& l_mix, Data& l_data) = 0;

};

#endif /* SETUPPROBLEM_H */
