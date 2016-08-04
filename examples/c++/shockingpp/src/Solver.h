#ifndef SOLVER_H
#define SOLVER_H

#include <boost/numeric/odeint.hpp>
#include <string>
#include "mutation++.h"

#include "SetupProperties.h"
#include "SetupProblem.h"

#include "Problem.h"
#include "Data.h"
#include "Mesh.h"
#include "BasicOutput.h"

using namespace boost::numeric::odeint;

class Solver {
public:
    Solver(size_t l_n_args, std::string l_file_name);
    ~Solver();

	void initialize();
	void solve();

private:
    size_t m_n_args;
    std::string s_file_name;
    std::string s_problem_type;
    SetupProperties* p_setup_props;
    SetupProblem* p_setup_problem;

    Mutation::Mixture* p_mix;
    Problem* p_problem;
    Mesh* p_mesh;
    ShockRelations* p_shock_relations;
    Data* p_data_pre;
    Data* p_data_post;
    BasicOutput* p_output;

    // This is for Larsen
    Data* p_data;

    inline void setMutationpp();
    inline void errorInsufficientNumberOfArguments();

};

#endif // SOLVER_H