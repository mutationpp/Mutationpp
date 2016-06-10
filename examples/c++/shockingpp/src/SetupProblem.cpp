#include "SetupProblem.h"
#include "SetupShocking1T.h"
#include "SetupShockingNT.h"

SetupProblem* SetupProblem::createFactory(std::string l_problem_type, std::string l_state_model){

    // @todo Replace this with switch and enum!
    if (!l_problem_type.compare("shocking") == true && !l_state_model.compare("ChemNonEq1T") == true) {
        return new SetupShocking1T();
    } else if(!l_problem_type.compare("shocking") == true && !l_state_model.compare("ChemNonEqTTv") == true) {
        return new SetupShockingNT();
    } else {
        std::cerr << "The requested type of problem has not been implemented yet!" << std::endl;
        exit(1);
    }

}


