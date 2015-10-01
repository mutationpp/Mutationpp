#include <iostream>
#include <stdlib.h>
#include <string>

#include "mutation++.h"

#include "InputData.h"
#include "InputFileParser.h"
#include "PostShockConditions.h"
#include "SolveSystem.h"

int main( int argc , char **argv ){

//=============================================================================================================

    // Check inputArguments
    std::string l_input_arg;
    const int default_test_case = 1;
    const std::string name_default_test_case( "shockair" );
    const int expected_number_arguments = 2;
    
    if ( argc == default_test_case ) {
        std::cout << "WARNING : No input file has been provided... The shockair.in input file will be loaded by default!" << std::endl;
        l_input_arg = name_default_test_case;
    } else if ( argc == expected_number_arguments ) {
        l_input_arg = argv[1];
    } else {
        std::cerr << "ERROR" << std::endl;
        exit(1);
    }

    // Opening input file and storing the data
    InputFileParser input_file_parser( l_input_arg );

    Mutation::Mixture* p_mixture = input_file_parser.getMixture();
    InputData* p_input_data = input_file_parser.getInputData();
    MeshOptions* p_mesh_options = input_file_parser.getMeshOptions();

    // Cross Shock Relations
    PostShockConditionsColdGas m_post_shock_conditions( p_mixture, p_input_data );

    // System
    SolveSystem m_solve_system( p_mixture, m_post_shock_conditions, p_mesh_options );

//    m_solve_system.setup();
//    m_solve_system.solve();

//=============================================================================================================

    return 0;

}
