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

  // Common for all the solvers
  static SetupProblem* createFactory(std::string l_problem_type, std::string l_state_model);
  virtual Problem* getProblem(Mutation::Mixture& l_mix, Data& l_data){};

  // S H O C K I N G 
  virtual Data* getDataPreShock(Mutation::Mixture& l_mix, const std::string& l_free_stream_conditions){};
  virtual Data* getDataPostShock(Mutation::Mixture& l_mix){};
  virtual ShockRelations* getShockRelations(Mutation::Mixture& l_mix){};

  // L A R S E N 
  virtual Data* getData(Mutation::Mixture& l_mix){};
};

#endif /* SETUPPROBLEM_H */
