#include "Solver.h"

int main(int argc, char** argv){

//=============================================================================================================

//@todo GIVE AN ELEGANT ERROR IN CASE ARGC == 0;
if(argc < 2) { std::cerr << "An input file should be provided! ABORTING!\n"; exit(1); }

    //@todo GIVE AN ELEGANT ERROR IN CASE ARGC == 0;
	Solver m_solver(argc, argv[1]);
	m_solver.initialize();
	m_solver.solve();

//=============================================================================================================

    return 0;
}
